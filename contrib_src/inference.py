import json
from processing import ImageProcessor
from modelhublib.model import ModelBase
from keras.models import model_from_json
from keras import backend as K
from keras.models import load_model
import tensorflow as tf

class Model(ModelBase):

    def __init__(self):
        # load config file
        config = json.load(open("model/config.json"))
        # get the image processor
        self._imageProcessor = ImageProcessor(config)
        # # load the DL model
        # load archiecture
        json_file = open( "model/architecture.json" , 'r')
        loaded_model_json = json_file.read()
        json_file.close()
        self._model = model_from_json(loaded_model_json)
        # load weights
        self._model.load_weights("model/weights.h5")
        self._graph = tf.get_default_graph()

    def infer(self, input):
        # load preprocessed input
        inputAsNpArr = self._imageProcessor.loadAndPreprocess(input)
        inputAsNpArrFormatted = inputAsNpArr.reshape(1,50,50,50,1)
        with self._graph.as_default():
            function = K.function([self._model.layers[0].input, K.learning_phase()],
            [self._model.layers[23].output])
            # Run inference with keras
            probabilities = self._model.predict_on_batch([inputAsNpArrFormatted])
            vector = function([inputAsNpArrFormatted, 0])
        # postprocess results into output
        output = self._imageProcessor.computeOutput([probabilities, vector[0]])
        return output
