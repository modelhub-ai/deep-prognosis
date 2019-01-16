from modelhublib.processor import ImageProcessorBase
import PIL
import SimpleITK
import numpy as np
import json


class ImageProcessor(ImageProcessorBase):

    def _preprocessBeforeConversionToNumpy(self, image):
        if isinstance(image, np.ndarray):
            image = image[50:100, 50:100, 50:100]
            image = self._centerAndNormalize(image)
        else:
            raise IOError("Image Type not supported for preprocessing.")
        return image

    def _preprocessAfterConversionToNumpy(self, npArr):
        # Since input is already a numpy array - there is nothing to do here.
        return npArr

    def computeOutput(self, inferenceResults):
        probs = np.squeeze(np.asarray(inferenceResults[0]))
        with open("model/labels.json") as jsonFile:
            labels = json.load(jsonFile)
        probLabel = []
        for i in range (len(probs)):
            obj = {'label': str(labels[str(i)]),
                    'probability': float(probs[i])}
            probLabel.append(obj)
        return [probLabel, inferenceResults[1].tolist()]

    def _centerAndNormalize(self, arr):
        out = arr
        oldMin = -1024
        oldRange = 3071+1024
        newRange = 1
        newMin = 0
        output = ((( out  - oldMin) * newRange * 1.0) / oldRange) + newMin
        return output
