import os
from modelhubapi import ModelHubAPI
from inference import Model

model = ModelHubAPI(Model(), os.path.dirname(os.path.realpath(__file__)))
