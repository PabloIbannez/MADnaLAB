import json

class MADnaModel:

    defaultModelData="./modelData/madna.json"

    def __init__(self,inputModelData:str=None):

        if inputModelData==None:
            self.modelData=self.defaultModelData
        else:
            self.modelData=self.inputModelData

        with open(self.modelData,"r") as f:
            self.model = json.load(f)

    def generateCoordinatesAndTopologyFromSeq(self,seq:str,output:str):
        return

if __name__=="__main__":
    m = MADnaModel()

