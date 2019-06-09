{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE TypeOperators #-}
{-# LANGUAGE FlexibleContexts #-}
{-# LANGUAGE ScopedTypeVariables #-}
{-# LANGUAGE ViewPatterns #-}

module Main where

import Control.Applicative (liftA2)
import Data.Array.Accelerate as A
import Data.Array.Accelerate.LLVM.Native as CPU
import Criterion.Main
import Control.DeepSeq


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
genVec :: Exp Int -> Acc (Vector Double)
genVec dim =
  generate (index1 dim) (const (A.fromIntegral $ 3 A.* dim - 1))

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
genMatrix' :: Exp Int -> Acc (Matrix Double)
genMatrix' dim = generate (index2 dim dim) f
  where
    f ((unlift . unindex2) -> (i,j)) =
      let fi :: Exp Double = A.fromIntegral i
       in A.cond (i A.== j)
            (fi A.+ A.fromIntegral dim)
            1.0

-- | Pack the LHS and RHS together so their dimensions always match
genMatrixVec' :: Exp Int -> Acc (Matrix Double, Vector Double)
genMatrixVec' = lift . liftA2 (,) genMatrix' genVec

genTest :: Exp Int -> Exp Int -> Acc (Vector Double, Scalar Double)
genTest k dim =
  let (m, ubvec) = unlift $ genMatrixVec' dim
      (r, d) = unlift . decompose $ m
      computation =
        aiterate (A.unit $ A.lift k) (jacobiIter r d ubvec) (fill (shape d) 0)
      error = error2 (mvm m computation) ubvec in
    A.lift (computation, error)

-- main :: IO ()
-- main = do
--   k :: Int <- readLn
--   -- read the number d, where the dimension dim = 2^d
--   dim :: Int <- (2 Prelude.^) <$> readLn
--   -- construct the computation and send everything to the CPU/GPU
--   print . CPU.run $ genTest 100 (2 Prelude.^ 14)

main :: IO ()
main
 = defaultMain
  [ bench "16 Iterations, 2^12 Dimensions" $ nf CPU.run (genTest 16 (2 Prelude.^ 12))
  , bench "100 Iterations, 2^12 Dimensions" $ nf CPU.run (genTest 100 (2 Prelude.^ 12))
  , bench "16 Iterations, 2^14 Dimensions" $ nf CPU.run (genTest 16 (2 Prelude.^ 14))
  ]

