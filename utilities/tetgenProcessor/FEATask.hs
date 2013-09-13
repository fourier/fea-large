module FEATask where

import Data.List
import Data.Tree

data Point = Point { xcoord :: Float
                   , ycoord :: Float
                   , zcoord :: Float
                   , nodeId :: Int
                   } deriving (Show)

type Element = [Int]

data SLAESolver = CG Float Int |
                  PCG_ILU Float Int |
                  Cholesky
    
type Attribute = (String, String)

data TaskLeaf = TreeLeaf String [Attribute] |
                   TreeNode Point |
                   TreeElement [Int]
                 deriving (Show)

type TaskTree = Tree TaskLeaf


stringFromTaskTree :: TaskTree -> String
stringFromTaskTree node =
    let root     = rootLabel node
        children = subForest node
        in
          "(" ++ foldl (\x y -> x ++ "\n" ++ y) (stringFromLeaf root) (map stringFromTaskTree children) ++ ")"
                        
    
stringFromLeaf :: TaskLeaf -> String
stringFromLeaf (TreeNode pnt) = pointAsSexp pnt
stringFromLeaf (TreeElement elt) = elementAsSexp elt
stringFromLeaf (TreeLeaf name attrs) =
    intercalate " " (name : map attributeAsString attrs)

readPoint :: String -> Point
readPoint str = let points = words str
                in
                  Point (read (points !! 1) :: Float)
                            (read (points !! 2) :: Float)
                            (read (points !! 3) :: Float)
                            (read (points !! 0) :: Int)

readElement :: String -> Element 
readElement str = map (\x -> (read x :: Int)-1) (drop 1 $ words str) 
                            
pointAsSexp :: Point -> String
pointAsSexp n = (show . xcoord $ n) ++ " " ++
                (show . ycoord $ n) ++ " " ++
                (show . zcoord $ n)

elementAsSexp :: Element -> String
elementAsSexp el = (intercalate " " (map show el))


-- returns the indicies of the nodes having specified x coordinate
ysectionPoints :: [Point] -> Float -> [Int]
ysectionPoints [] _ = []
ysectionPoints points yc =
    map nodeId $ filter (\(Point _ y _ _) -> y == yc) points
        
-- returns indices of the nodes by criterion
leftmostPoints :: [Point] -> [Int]
leftmostPoints points = 
    let (Point _ yc _ _) = minimumBy (\a b -> compare (ycoord a)
                                              (ycoord b)) points
    in
      ysectionPoints points yc

-- returns indices of the rightmost nodes
rightmostPoints :: [Point] -> [Int]
rightmostPoints points = 
    let (Point _ yc _ _) = maximumBy (\a b -> compare (ycoord a)
                                              (ycoord b)) points
    in
      ysectionPoints points yc
          
                   
attribute :: String -> String -> Attribute
attribute name value = (name,value)
                     
attributeAsString :: Attribute -> String
attributeAsString (name,value) = ":" ++ name ++ " " ++ value

nodeWithAttributes :: String -> [Attribute] -> TaskTree
nodeWithAttributes name attrs =
    let
        attrsList = map (\(x,y) -> attribute x y) attrs
    in
      Node (TreeLeaf name attrsList) []

                                
modelParameters :: Int -> Int -> TaskTree
modelParameters mu lambda = nodeWithAttributes "model-parameters"
                            [("mu", show mu),
                             ("lambda", show lambda)]

    
lineSearch :: Int -> TaskTree
lineSearch count = nodeWithAttributes "line-search" [("max", show count)]

arcLength :: Int -> TaskTree
arcLength count = nodeWithAttributes "arc-length" [("max", show count)]

                
slaeSolver :: SLAESolver -> TaskTree
slaeSolver Cholesky = nodeWithAttributes "slae-solver" [("type","CHOLESKY")]
slaeSolver (CG tol iter) =
    nodeWithAttributes "slae-solver" [("type","CG"),
                                      ("tolerance",show tol),
                                      ("max-iterations",show iter)]
slaeSolver (PCG_ILU tol iter) =
    nodeWithAttributes "slae-solver" [("type","PCG_ILU"),
                                      ("tolerance",show tol),
                                      ("max-iterations",show iter)]
                  
elementType :: String -> Int -> Int -> TaskTree
elementType name gaussCount nodesCount =
    nodeWithAttributes "element-type"
                       [("gauss-nodes-count", show gaussCount),
                        ("name", name),
                        ("nodes-count", show nodesCount)]

model :: String -> TaskTree
model name = let mdl = nodeWithAttributes "model" [("name", name)]
             in Node (rootLabel mdl) [modelParameters 100 100]

solution :: Float -> Int -> Int -> Bool -> TaskTree
solution tolerance increments maxNewton modified =
    let mdl = nodeWithAttributes "solution"
              [("desired-tolerance", show tolerance),
               ("task-type", "CARTESIAN3D"),
               ("load-increments-count", show increments),
               ("modified-newton", (if (modified)
                                    then "yes"
                                    else "no")),
               ("max-newton-count", show maxNewton)]
    in
      Node (rootLabel mdl) [elementType "TETRAHEDRA10" 5 10,
                            slaeSolver (CG 1e-14 15000),
                            lineSearch 0,
                            arcLength 0]


--solution 1e-6 10 100 True 

nodesTree :: [Point] -> TaskTree
nodesTree arr = Node (TreeLeaf "nodes" [])
                (map (\x -> Node (TreeNode x) []) arr)


elementsTree :: [Element] -> TaskTree
elementsTree arr = Node (TreeLeaf "elements" [])
                   (map (\x -> Node (TreeElement x) []) arr)

geometry :: [Point] -> [Element] -> TaskTree
geometry nodes elements = Node (TreeLeaf "geometry" [])
                          [nodesTree nodes,
                           elementsTree elements]


prescNode :: Float -> Float -> Float -> Int -> Int -> TaskTree
prescNode xoffset yoffset zoffset bndType index = 
    nodeWithAttributes "presc-node"
                       [("x", show xoffset),
                        ("y", show yoffset),
                        ("z", show zoffset),
                        ("type", show bndType),
                        ("node-id", show index)]

prescribedDisplacements :: [Int] -> [Int] -> Float -> TaskTree
prescribedDisplacements (left:leftmost) rightmost offset =
                           let mdl = nodeWithAttributes "prescribed-displacements" []
                               leftNodes = map (nodeWithOffset 0) leftmost
                               nodes = (prescNode 0 0 0 7 left:leftNodes) ++ (map (nodeWithOffset 0.05) rightmost)
                           in
                             Node (rootLabel mdl)
                                  nodes
                           where 
                             nodeWithOffset offset idx =
                                 (prescNode 0 offset 0 2 idx)

boundaryConditions :: [Int] -> [Int] -> TaskTree
boundaryConditions leftmost rightmost =
    Node (TreeLeaf "boundary-conditions" [])
             [prescribedDisplacements leftmost rightmost 0.05]


--- gathering all together
finiteElementTaskTree :: [Point] -> [Element] -> TaskTree
finiteElementTaskTree nodes elements =
    let leftmostIndices = leftmostPoints nodes
        rightmostIndices = rightmostPoints nodes
        modelTree = model "COMPRESSIBLE_NEOHOOKEAN"
        solutionTree = solution 1e-6 100 100 True
        inputDataTree = Node (TreeLeaf "input-data" [])
                        [geometry nodes elements, boundaryConditions leftmostIndices rightmostIndices]
    in
      Node (TreeLeaf "task" [])
           [modelTree,solutionTree,inputDataTree]
                              
