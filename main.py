# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 13:34:11 2022

@author: maria
"""
from neuron import h
h("lsd = 0") # 0: no lsd, 1: on a giant trip, anything in between 0 and 1 = smaller trip
h.load_file("stdlib.hoc")
h.load_file("nrngui.hoc")
h.load_file("altered_complexity_model.hoc")
