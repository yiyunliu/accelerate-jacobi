{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE ApplicativeDo #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE TemplateHaskell #-}
module Main where
import Data.Array.Accelerate as A
import Data.Array.Accelerate.LLVM.Native as A
import Control.Applicative (liftA2)

genMatrix :: Int -> Matrix Double
genMatrix dim = fromFunction (Z :. dim :. dim) f
  where f (_ :. i :. j) =
          let fi :: Double = Prelude.fromIntegral i in
          if i Prelude.== j
          then 0 -- fi + 1.0
          else 1.0

genDiagonal :: Int -> Vector Double
genDiagonal dim = fromFunction (Z :. dim) f
  where f (_ :. i) =
          let fi :: Double = Prelude.fromIntegral i in fi + 1.0
            
genVec :: Int -> Vector Double
genVec dim = fromFunction (Z :. dim)
  (const (Prelude.fromIntegral $ 3 Prelude.* dim - 1))


genMatrixVec :: Int -> (Matrix Double, Vector Double, Vector Double)
genMatrixVec = do
  m <- genMatrix
  d <- genDiagonal
  v <- genVec
  pure (m, d, v)

dotp :: A.Num a => Acc (Vector a) -> Acc (Vector a) -> Acc (Scalar a)
dotp v  = A.fold (A.+) 0 . A.zipWith (*) v

error2 :: A.Floating a => Acc (Vector a) -> Acc (Vector a) -> Acc (Scalar a)
error2 v = A.map A.sqrt . A.fold (A.+) 0 . A.map (A.** 2) . A.zipWith (-) v

mvm :: A.Num a => Acc (Matrix a) -> Acc (Vector a) -> Acc (Vector a)
mvm mat vec =
  let Z :. rows :. cols = unlift (shape mat) :: Z :. Exp Int :. Exp Int
      vec'              = A.replicate (lift (Z :. rows :. All)) vec
  in
  A.fold (+) 0 ( A.zipWith (*) mat vec' )


jacobiIter :: Acc (Matrix Double) -> Acc (Vector Double) ->
  Acc (Vector Double) -> Acc (Vector Double) -> Acc (Vector Double)
jacobiIter m d v x = A.zipWith (/) (A.zipWith (A.-) v (mvm m x)) d


aiterate
  :: forall a. Arrays a -- forall is necessary
  => Acc (Scalar Int)
  -> (Acc a -> Acc a)
  -> Acc a
  -> Acc a
-- scoped type variable
aiterate n f z
  = let step :: (Acc (Scalar Int, a)) -> (Acc (Scalar Int, a))
        step (unlift -> (i, acc))   = lift (A.map succ i, f acc)
    in
    A.asnd $ awhile (\v -> A.zipWith (A.<) (A.afst v)  n) step (lift (unit 0, z))


matrix :: Matrix Double
matrix = fromList (Z :. 3 :. 3) [5, -2, 3, -3, 9, 1, 2, -1, -7]

bvec :: Vector Double
bvec = fromList (Z :. 3) [-1, 2 ,3]

decompose :: Acc (Matrix Double) -> Acc (Matrix Double, Vector Double)
decompose m = lift (imap g m, d)
  where d = let dim = indexTail (shape m) in
          generate dim f
        f :: Exp (Z :. Int) -> Exp Double
        f (unindex1 -> i) = m ! (index2 i i)
        g :: Exp (Z :. Int :. Int) -> Exp Double -> Exp Double
        g ((unlift . unindex2) -> (i, j)) e = A.cond (i A.== j) 0 e


genMatrix' :: Int -> Matrix Double
genMatrix' dim = fromFunction (Z :. dim :. dim) f
  where f (_ :. i :. j) =
          let fi :: Double = Prelude.fromIntegral i in
          if i Prelude.== j
          then fi Prelude.+ Prelude.fromIntegral dim
          else 1.0

genMatrixVec' :: Int -> (Matrix Double, Vector Double)
genMatrixVec' = liftA2 (,) genMatrix' genVec


-- main :: IO ()
-- main = do
--   k :: Int <- readLn
--   dim :: Int <- (2 Prelude.^) <$> readLn
--   let (m, d, v) = genMatrixVec dim
--   let computation = aiterate (A.unit $ A.lift k) (jacobiIter (use m) (use d) (use v) )  (use $ fromList (Z :. dim) (repeat 0))
--   print . run $ computation
  

-- main :: IO ()
-- main = do
--   k :: Int <- readLn
--   let m = use matrix
--       (r,d) = unlift . decompose $ m
--       ubvec = use bvec
--       computation = aiterate (A.unit $ A.lift k) (jacobiIter r d ubvec)  (fill (shape d) 0)
--       error = error2 (mvm m computation) ubvec
--   print $ A.run (lift (computation, error))

main :: IO ()
main = do
  k :: Int <- readLn
  dim :: Int <- (2 Prelude.^) <$> readLn
  let (matrix, bvec) = genMatrixVec' dim
  -- print matrix
  let m = use matrix
      (r,d) = unlift . decompose $ m
      ubvec = use bvec
      computation = aiterate (A.unit $ A.lift k) (jacobiIter r d ubvec)  (fill (shape d) 0)
      error = error2 (mvm m computation) ubvec
  print $ A.run (lift (computation, error))
