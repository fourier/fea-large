module Main( main ) where

import System.Environment( getArgs )
import System.FilePath.Posix
import FEATask
    
                    
                    
namesFromFilename :: String -> (String,String,String)
namesFromFilename f =
    let (fname,_) = splitExtension f
    in
      (fname ++ ".sexp", fname ++ ".1.node", fname ++ ".1.ele")

skipSpaces = dropWhile (\x -> (x == ' ') || ( x == '\t'))
      
      
filterLines =
    filter (\x -> (notEmpty x) && (notComment x))
    where
      notEmpty x = (length $ words x) /= 0
      notComment x = let contents = words x
                     in
                       case contents of
                         []     -> False
                         (x:xs) -> (head x) /= '#'

nodesVerifyLines filteredLines =
    (read . head . words . head $ filteredLines :: Int) == (length filteredLines) - 1


readDataLines :: (String -> a) -> [String] -> [a]
readDataLines readFunc lines  = 
    if (nodesVerifyLines lines) then
        map readFunc (drop 1 lines)
    else
        error "Incorrect number of lines in file"
                                                        
                                                        
readPoints :: [String] -> [Point]
readPoints = readDataLines readPoint 

readElements :: [String] -> [Element]
readElements = readDataLines readElement
              
readFiles nodename elename =
    do
      nodeContents <- readFile nodename
      eleContents  <- readFile elename
      return (filterLines . lines $ nodeContents,
              filterLines . lines $ eleContents)


processFiles outname nodename elename =
    do
      (nodeLines,eleLines) <- readFiles nodename elename
      let nodes = readPoints nodeLines
      let elements = readElements eleLines
      putStrLn ("Number of nodes: " ++ show (length nodes ))
      putStrLn ("Number of elements: " ++ show (length elements))
      writeFile outname (";; -*- Mode: lisp; -*-\n" ++ (stringFromTaskTree $ finiteElementTaskTree nodes elements))
      putStrLn "Done!"

main = do
  args <- getArgs
  case args of
    []  -> error "Filename (like brick.poly) needed"
    [f] -> let (outname,nodename,elename) = namesFromFilename f
           in
             processFiles outname nodename elename
                     
