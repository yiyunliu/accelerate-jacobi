{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE ApplicativeDo #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ViewPatterns #-}
{-# LANGUAGE TemplateHaskell #-}

module Main where

import Control.Applicative (liftA2)
import Data.Array.Accelerate as A
import Data.Array.Accelerate.LLVM.Native as A

-- | Dot product. stolen from accelerate's doc
dotp :: A.Num a => Acc (Vector a) -> Acc (Vector a) -> Acc (Scalar a)
dotp v = A.fold (A.+) 0 . A.zipWith (*) v

-- | Calculating error by 2-norm
error2 :: A.Floating a => Acc (Vector a) -> Acc (Vector a) -> Acc (Scalar a)
error2 v = A.map A.sqrt . A.fold (A.+) 0 . A.map (A.** 2) . A.zipWith (-) v

-- | Matrix multiplication. stolen from accelerate's doc
mvm :: A.Num a => Acc (Matrix a) -> Acc (Vector a) -> Acc (Vector a)
mvm mat vec =
  let Z :. rows :. cols = unlift (shape mat) :: Z :. Exp Int :. Exp Int
      vec' = A.replicate (lift (Z :. rows :. All)) vec
   in A.fold (+) 0 (A.zipWith (*) mat vec')

-- | One iteration of jacobi's method
jacobiIter ::
     Acc (Matrix Double)
  -> Acc (Vector Double)
  -> Acc (Vector Double)
  -> Acc (Vector Double)
  -> Acc (Vector Double)
jacobiIter m d v x = A.zipWith (/) (A.zipWith (A.-) v (mvm m x)) d

-- | Generate a hardwired RHS of the linear equation
genVec :: Int -> Vector Double
genVec dim =
  fromFunction (Z :. dim) (const (Prelude.fromIntegral $ 3 Prelude.* dim - 1))

-- | Acc version of iterate
aiterate ::
     forall a. Arrays a -- forall is necessary
  => Acc (Scalar Int)
  -> (Acc a -> Acc a)
  -> Acc a
  -> Acc a
aiterate n f z =
  let step :: (Acc (Scalar Int, a)) -> (Acc (Scalar Int, a)) -- scoped type variable
      step (unlift -> (i, acc)) = lift (A.map succ i, f acc)
   in A.asnd $
      awhile (\v -> A.zipWith (A.<) (A.afst v) n) step (lift (unit 0, z))

-- | Decompose a matrix M into D (represented as vector) and L + U
decompose :: Acc (Matrix Double) -> Acc (Matrix Double, Vector Double)
decompose m = lift (imap g m, d)
  where
    d =
      let dim = indexTail (shape m)
       in generate dim f
    f :: Exp (Z :. Int) -> Exp Double
    f (unindex1 -> i) = m ! (index2 i i)
    g :: Exp (Z :. Int :. Int) -> Exp Double -> Exp Double
    g ((unlift . unindex2) -> (i, j)) e = A.cond (i A.== j) 0 e

-- | Generate a hardwired diagonally dominant LHS of the linear equation.
genMatrix' :: Int -> Matrix Double
genMatrix' dim = fromFunction (Z :. dim :. dim) f
  where
    f (_ :. i :. j) =
      let fi :: Double = Prelude.fromIntegral i
       in if i Prelude.== j
            then fi Prelude.+ Prelude.fromIntegral dim
            else 1.0

-- | Pack the LHS and RHS together so their dimensions always match
genMatrixVec' :: Int -> (Matrix Double, Vector Double)
genMatrixVec' = liftA2 (,) genMatrix' genVec

main :: IO ()
main
  -- read the number of iterations
 = do
  k :: Int <- readLn
  -- read the number d, where the dimension dim = 2^d
  dim :: Int <- (2 Prelude.^) <$> readLn
  let (matrix, bvec) = genMatrixVec' dim
  -- construct the computation before sending everything to the CPU/GPU
  let m = use matrix
      (r, d) = unlift . decompose $ m
      ubvec = use bvec
      computation =
        aiterate (A.unit $ A.lift k) (jacobiIter r d ubvec) (fill (shape d) 0)
      error = error2 (mvm m computation) ubvec
  -- run the program and print out the result as a tuple containing the solution and error
  print $ A.run (lift (computation, error))
