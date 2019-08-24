function pushDownLocal(atTree::BallTreeDensity, aRoot::Int, hdl::pArrHdls)
    if !(isLeaf(atTree, aRoot))
      close = atTree.left(aRoot);
      if (close != NO_CHILD)
        hdl.pAdd[1,close] += hdl.pAdd[1,aRoot]
      end
      close = right(atTree, aRoot);
      if (close != NO_CHILD)
        hdl.pAdd[1,close] += hdl.pAdd[1,aRoot]
      end
      hdl.pAdd[1,aRoot] = 0
    end
end

function pushDownAll(locations::BallTreeDensity, hdl::pArrHdls)
  for j in root():(leafFirst(locations,root())-1)
      hdl.pAdd[1,  left(locations, j)] += hdl.pAdd[1,j]
      hdl.pAdd[1, right(locations, j)] += hdl.pAdd[1,j]
      hdl.pAdd[1,j] = 0
    end
    for j in leafFirst(locations, root()):(leafLast(locations, root())+1)
      hdl.pMin[j] += hdl.pAdd[1,j] - hdl.pErr[j]
      hdl.pMax[j] += hdl.pAdd[1,j] + hdl.pErr[j]
      hdl.pAdd[1,j] = 0; hdl.pErr[j] = 0;
    end
end
