class Node():
    def __init__(self):
        self.parent = None
        self.leaf = False
        self.left_child = None
        self.right_child = None
        self.overlap = None
        self.center = None
        self.interval_info = None

class Interval_tree():
    """docstring for Interval_tree, only use integers, listed intervals should be in [start,end,interval_info] format bounds inclusive"""
    def __init__(self):
        self.root = None
        self.list = None
        self.startList = None
        self.stopList = None
        pass
    def build_from_list(self,list,start,end,parent=None):
        this_node = Node()
        if parent==None:
            self.list = list
            self.root = this_node
        this_node.center = (end-start)//2 + start
        this_node.parent = parent
        s_left = []
        s_right = []
        this_node.overlap = []
        for interval in list:
            if interval[1] < this_node.center:
                s_left.append(interval)
            elif interval[0] > this_node.center:
                s_right.append(interval)
            else:
                this_node.overlap.append(interval)
        if len(s_left)>0:
            this_node.left_child = self.build_from_list(s_left,start,this_node.center,parent=this_node)
        if len(s_right)>0:
            this_node.right_child = self.build_from_list(s_right,this_node.center,end,parent=this_node)
        if parent==None:
            self.startBST = self.buildBST(sorted(self.list,key=lambda x:x[0]))
            self.endBST = self.buildBST(sorted(self.list,key=lambda x:x[1]))
        return this_node

    def buildBST(self,list,start=None,stop=None,parent=None):
        if parent==None:
            start = 0
            stop = len(list)
        diff = stop-start
        this_node = Node()
        this_node.parent = parent
        if(diff>1):
            mid = diff//2 + start
            this_node.interval_info = list[mid]
            this_node.left_child = self.buildBST(list,start,mid,this_node)
            this_node.right_child = self.buildBST(list,mid+1,stop,this_node)
        elif diff>0:
            this_node.interval_info = list[start]

    def query(self,val,node="root"):
        if node == "root": node = self.root
        if not node: return []
        intersecting = None
        if val==node.center:
            return node.overlap[:]
        elif(val<node.center):
            intersecting = [interval for interval in node.overlap if interval[0]<val]
            intersecting+= self.query(val,node.left_child)
        else:
            intersecting = [interval for interval in node.overlap if interval[1]>val]
            intersecting+= self.query(val,node.right_child)
        return intersecting

    def lines(self,node="root"):
        if node == "root": node = self.root
        if not node: return []
        linesList = [str(node.center)+": "]#+str(node.overlap)]
        linesList+= ["\t"+line for line in self.lines(node.left_child)]
        linesList+= ["\t"+line for line in self.lines(node.right_child)]
        return linesList