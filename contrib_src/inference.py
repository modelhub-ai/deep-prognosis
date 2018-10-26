import onnx
import caffe2.python.onnx.backend
import json
from processing import ImageProcessor
from modelhublib.model import ModelBase


class Model(ModelBase):

    def __init__(self):
        # load config file
        config = json.load(open("model/config.json"))
        # get the image processor
        self._imageProcessor = ImageProcessor(config)
        # load the DL model (change this if you are not using ONNX)
        self._model = onnx.load('model/model.onnx')
    

    def infer(self, input):
        # load preprocessed input
        inputAsNpArr = self._imageProcessor.loadAndPreprocess(input)
        # Run inference with caffe2 (change this if you are using a different DL framework)
        results = caffe2.python.onnx.backend.run_model(self._model, [inputAsNpArr])
        # postprocess results into output
        output = self._imageProcessor.computeOutput(results)
        return output
        

