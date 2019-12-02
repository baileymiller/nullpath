# Author: Bailey Miller
# Description: Run parallel pbrt jobs on anthill cluster.

import argparse
import ast
import json
import os
import subprocess
import time
import itertools
import math
from copy import deepcopy

class FLAG_TYPES:
  INTEGRATOR = 'integrator'

class COMMON_FLAGS:
  enable_equal_majorant = {
    "type": "bool",
    "name": "enable-equal-majorant",
    "value": "\"true\""
  }
  enable_power_heuristic = {
    "type": "bool",
    "name": "enable-balance-heuristic",
    "value": "\"false\""
  }
  enable_nee = {
    "type": "bool", 
    "display": "strategy=dir-mis",
    "name": "enable-NEE",
    "value": "\"true\"",
  }
  enable_only_nee = {
    "type": "bool", 
    "display": "strategy=nee",
    "name": "enable-only-NEE",
    "value": "\"true\"",
  }
  enable_only_ea = {
    "type": "bool", 
    "display": "strategy=only_ea",
    "name": "enable-only-EA",
    "value": "\"true\"",
  }
  enable_dist_mis = {
    "type": "bool", 
    "display": "strategy=dist-mis",
    "name": "enable-EA",
    "value": "\"true\"",
  }
  enable_split_nee = {
    "type": "bool",
    "display": "strategy=mis-nee",
    "name": "enable-split-NEE",
    "value": "\"true\"",
  }
  enable_light_sampling = {
    "type": "bool",
    "display": "strategy=light-sampler",
    "name": "enable-LightSample",
    "value": "\"true\"",
  } 
  enable_naive = { 
    "display": "strategy=naive",
  }
  enable_naive_delta = { 
    "type": "bool",
    "display": "strategy=naive-delta",
    "name": "enable-naive-delta",
    "value": "\"true\""
  }
  enable_delta = { 
    "type": "bool",
    "display": "strategy=delta",
    "name": "enable-delta-track",
    "value": "\"true\""
  }
  single_scatter = {
    "type": "bool", 
    "display": "ss",
    "name": "singleScatter",
    "value": "\"true\""
  }
  ignore_no_scatter = {
    "type": "bool",
    "display": "at-least-one-scatter",
    "name": "ignoreNoScatter",
    "value": "\"true\""
  }
  print_path = {
    "type": "bool",
    "display": None,
    "name": "printPath",
    "value": "\"true\""
  }
  print_hits = {
    "type": "bool",
    "display": None,
    "name": "printHits",
    "value": "\"true\""
  }
  print_ea_path = {
    "type": "bool",
    "display": None,
    "name": "printEAPath",
    "value": "\"true\""
  }
  show_nee_mis_weights = {
    "type": "bool",
    "display": "mis-weights",
    "name": "showTwoStrategyMISWeight",
    "value": "\"true\""
  }
  split_nee_mis_values = {
    "type": "bool",
    "display": "mis-values",
    "name": "splitTwoStrategyMISValue",
    "value": "\"true\""
  }
  show_dir_pdf = {
    "type": "bool",
    "display": "dirPdf",
    "name": "showDirPdf",
    "value": "\"true\""
  }
  enable_surface_nee = {
    "type": "bool",
    "display": "surfaceNEE",
    "name": "enable-surface-NEE",
    "value": "\"true\""
  }
  @staticmethod
  def maxDepth(x):
    return {
      "type": "integer",
      "display": "maxDepth={0}".format(x),
      "name": "maxdepth",
      "value": "[{0}]".format(x)
    }

  @staticmethod
  def constantEmission(x): 
    return {
      "type": "float", 
      "display": "vol-emission={0}".format(x),
      "name": "constantEmission",
      "value": "[{0}]".format(x)
    }
  @staticmethod
  def viewNumScatterEvents(x): 
    return {
      "type": "integer", 
      "display": "num-scatter-events={0}".format(x),
      "name": "viewNumScatterEvents",
      "value": "[{0}]".format(x)
    }

  @staticmethod
  def set_p_channel(w_r, w_g, w_b):
    w_total = float(w_r + w_g + w_b)
    return {
      "type": "rgb",
      "display": "channel-probability=[{0}, {1}, {2}]".format(w_r / w_total, 
                                        w_g / w_total, 
                                        w_b / w_total),
      "name": "p_channel",
      "value": "[{0} {1} {2}]".format(w_r / w_total, 
                                        w_g / w_total, 
                                        w_b / w_total)
    }

# EXPERIMENTS
def setup_veach_comparison_one(args):
  args.xres = 800
  args.yres = 800
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure3/scene_template.pbrt"]
  args.integrators = ["nullpath"] #, "nullpath"]
  args.integrator_flags = {
    "nullpath" : [
        {"flags" : [COMMON_FLAGS.enable_naive,
                    COMMON_FLAGS.single_scatter,
                    COMMON_FLAGS.ignore_no_scatter],
         "samples" : 1},
        #{"flags" : [COMMON_FLAGS.enable_naive,
        #            COMMON_FLAGS.single_scatter,
        #            COMMON_FLAGS.ignore_no_scatter],
        #  "samples" : 32},
        #{"flags" : [COMMON_FLAGS.enable_only_nee,
        #            COMMON_FLAGS.single_scatter,
        #            COMMON_FLAGS.ignore_no_scatter],
        #  "samples" : 30},
        #{"flags" : [COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.single_scatter,
        #            COMMON_FLAGS.ignore_no_scatter],
        #  "samples": 27}
      ],
      "spectral" : [
         #{"flags" : [COMMON_FLAGS.enable_naive,
         #            COMMON_FLAGS.single_scatter,
         #            COMMON_FLAGS.ignore_no_scatter],
         # "samples" : 32},
         #{"flags" : [COMMON_FLAGS.enable_only_nee,
         #            COMMON_FLAGS.single_scatter,
         #            COMMON_FLAGS.ignore_no_scatter],
         # "samples" : 31},
         #{"flags" : [COMMON_FLAGS.enable_split_nee,
         #            COMMON_FLAGS.single_scatter,
         #            COMMON_FLAGS.ignore_no_scatter],
         #    "samples": 28},
        ]
  }
  args.custom_scene_file = True
  args.custom_formatting = []
  
  naive_opts = {
    'scene_id' : 'veach_one', 
    'g' : 0.99, 
    'r' : 3,
    'translate': '-4.6 0 0',
    's' : 4,
    'a' : 0,
    'n' : 0,
    'l' : 1,
    'gamma': 1,
  }
  args.custom_formatting.append(naive_opts)
  args.threads = 8
  return args

def setup_veach_comparison_two(args):
  args.xres = 800
  args.yres = 800
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure3/scene_template.pbrt"]
  args.integrators = ["nullpath"] #, "spectral"]
  args.integrator_flags = {
      "nullpath" : [
        #{"flags" : [COMMON_FLAGS.enable_only_nee,
        #            COMMON_FLAGS.single_scatter,
        #            COMMON_FLAGS.ignore_no_scatter],
        # "samples" : 3000},
        {"flags" : [COMMON_FLAGS.enable_naive,
                    COMMON_FLAGS.single_scatter,
                    COMMON_FLAGS.ignore_no_scatter],
         "samples" : 32},
        #{"flags" : [COMMON_FLAGS.enable_only_nee,
         #           COMMON_FLAGS.single_scatter,
         #           COMMON_FLAGS.ignore_no_scatter],
         #"samples" : 26},
        #{"flags" : [COMMON_FLAGS.enable_split_nee,
        #             COMMON_FLAGS.single_scatter,
        #             COMMON_FLAGS.ignore_no_scatter],
        #  "samples": 23},
      ],
      "spectral" : [
          #{"flags" : [COMMON_FLAGS.enable_naive,
          #           COMMON_FLAGS.single_scatter,
          #           COMMON_FLAGS.ignore_no_scatter],
          #  "samples" : 32},
          #{"flags" : [COMMON_FLAGS.enable_only_nee,
          #            COMMON_FLAGS.single_scatter,
          #            COMMON_FLAGS.ignore_no_scatter],
          #  "samples": 32},
          {"flags" : [COMMON_FLAGS.enable_split_nee,
                      COMMON_FLAGS.single_scatter,
                      COMMON_FLAGS.ignore_no_scatter],
             "samples": 28},
        ]
  }

  args.custom_scene_file = True
  args.custom_formatting = []
  
  nee_opts = {
    'scene_id' : 'veach_two_ss', 
    'g' : 0.99, 
    'r' : 0.3,
    'translate': '0 2 2',
    's' : 16,
    'a' : 0,
    'n' : 0,
    'gamma': 1,
    'l' : 10000,
  }
  args.custom_formatting.append(nee_opts)
  args.threads = 8
  return args

def setup_colored_smoke_one(args):
  args.xres = 800
  args.yres = 1600
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure5/smoke_template.pbrt"]
  args.integrators = ["spectral"]
  args.integrator_flags = {
      "spectral" : [
          {"flags": [COMMON_FLAGS.enable_split_nee,
                     COMMON_FLAGS.single_scatter,
                     COMMON_FLAGS.ignore_no_scatter],
           "samples": 3000},
          #{"flags": [COMMON_FLAGS.enable_split_nee,
          #           COMMON_FLAGS.enable_naive_delta,
          #           COMMON_FLAGS.single_scatter,
          #           COMMON_FLAGS.ignore_no_scatter],
          #           "samples": 33},
          #{"flags": [COMMON_FLAGS.enable_split_nee,
          #           COMMON_FLAGS.single_scatter,
          #           COMMON_FLAGS.ignore_no_scatter],
          #           "samples": 17},
       ],
      "nullpath" : [
        {"flags" : [COMMON_FLAGS.enable_delta,
                    COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.single_scatter,
                    COMMON_FLAGS.ignore_no_scatter],
            "samples" : 32},
        #{"flags" : [COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.single_scatter,
        #            COMMON_FLAGS.ignore_no_scatter],
        # "samples" : 32},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  noise_x = 0.5
  noise_y = 4
  noise_z = 0.5

  density_r = 30
  albedo_r = 0.001

  density_g = 100
  albedo_g = 0.001

  density_b = 15
  albedo_b = 0.95

  density_r_2 = 15
  albedo_r_2 = 0.95

  density_g_2 = 100
  albedo_g_2 = 0.001

  density_b_2 = 30
  albedo_b_2 = 0.001

  opts = {
    'scene_id': 'smoke_one_ss', 
    'g' : -0.4, 
    'l' : 0.75,
    'gamma': 0.15,
    's_r' : density_r * albedo_r, 
    's_g' : density_g * albedo_g, 
    's_b' : density_b * albedo_b, 
    'a_r' : density_r * (1.0 - albedo_r),
    'a_g' : density_g * (1.0 - albedo_g),
    'a_b' : density_b * (1.0 - albedo_b),
    's_r_2' : density_r_2 * albedo_r_2, 
    's_g_2' : density_g_2 * albedo_g_2, 
    's_b_2' : density_b_2 * albedo_b_2, 
    'a_r_2' : density_r_2 * (1.0 - albedo_r_2), 
    'a_g_2' : density_g_2 * (1.0 - albedo_g_2),
    'a_b_2' : density_b_2 * (1.0 - albedo_b_2),
    'noise_x': noise_x,
    'noise_y': noise_y,
    'noise_z': noise_z,
  }

  args.custom_formatting.append(opts)
  args.threads = 16
  return args

def setup_colored_smoke_one_ms_new(args):
  args.xres = 800
  args.yres = 1600
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure5/smoke_template.pbrt"]
  args.integrators = ["nullpath"] #, "nullpath"]
  args.integrator_flags = {
      "spectral" : [
          #{"flags": [COMMON_FLAGS.enable_split_nee,
          #           COMMON_FLAGS.ignore_no_scatter],
          # "samples": 3000},
          {"flags": [COMMON_FLAGS.enable_split_nee,
                     COMMON_FLAGS.enable_naive_delta,
                     COMMON_FLAGS.ignore_no_scatter],
                     "samples": 32},
         #{"flags": [COMMON_FLAGS.enable_split_nee,
         #            COMMON_FLAGS.ignore_no_scatter],
         #            "samples": 16},
       ],
      "nullpath" : [
        {"flags" : [COMMON_FLAGS.enable_delta,
                    COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.ignore_no_scatter],
            "samples" : 30},
        #{"flags" : [COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.ignore_no_scatter],
        #    "samples" : 27},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  noise_x = 0.5
  noise_y = 4
  noise_z = 0.5

  density_r = 30
  albedo_r = 0.95

  density_g = 100
  albedo_g = 0.001

  density_b = 15
  albedo_b = 0.001

  density_r_2 = 15
  albedo_r_2 = 0.001

  density_g_2 = 100
  albedo_g_2 = 0.001

  density_b_2 = 30
  albedo_b_2 = 0.95

  opts = {
    'scene_id': 'smoke_one_ss', 
    'g' : -0.4, 
    'l' : 0.75,
    'gamma': 0.15,
    's_r' : density_r * albedo_r, 
    's_g' : density_g * albedo_g, 
    's_b' : density_b * albedo_b, 
    'a_r' : density_r * (1.0 - albedo_r),
    'a_g' : density_g * (1.0 - albedo_g),
    'a_b' : density_b * (1.0 - albedo_b),
    's_r_2' : density_r_2 * albedo_r_2, 
    's_g_2' : density_g_2 * albedo_g_2, 
    's_b_2' : density_b_2 * albedo_b_2, 
    'a_r_2' : density_r_2 * (1.0 - albedo_r_2), 
    'a_g_2' : density_g_2 * (1.0 - albedo_g_2),
    'a_b_2' : density_b_2 * (1.0 - albedo_b_2),
    'noise_x': noise_x,
    'noise_y': noise_y,
    'noise_z': noise_z,
  }

  args.custom_formatting.append(opts)
  args.threads = 8
  return args


def setup_colored_smoke_two(args):
  args.xres = 800
  args.yres = 1600
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure5/smoke_template.pbrt"]
  args.integrators = ["nullpath"] #, "nullpath"]
  args.integrator_flags = {
      "spectral" : [
         #{"flags": [COMMON_FLAGS.enable_split_nee,
         #           COMMON_FLAGS.single_scatter,
         #           COMMON_FLAGS.ignore_no_scatter],
         #            "samples": 3000},
         {"flags": [COMMON_FLAGS.enable_naive_delta,
                    COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.single_scatter,
                    COMMON_FLAGS.ignore_no_scatter],
                     "samples": 35},
         # {"flags": [COMMON_FLAGS.enable_split_nee,
         #            COMMON_FLAGS.single_scatter,
         #            COMMON_FLAGS.ignore_no_scatter],
         #            "samples": 20},
       ],
      "nullpath" : [
        #{"flags" : [COMMON_FLAGS.enable_delta,
        #            COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.single_scatter,
        #            COMMON_FLAGS.ignore_no_scatter],
        #    "samples" : 35},
        {"flags" : [COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.single_scatter,
                 COMMON_FLAGS.ignore_no_scatter],
          "samples" : 32},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  noise_x = 0.5
  noise_y = 2
  noise_z = 0.5

  density_r = 30
  albedo_r = 0.001

  density_g = 100
  albedo_g = 0.95

  density_b = 15
  albedo_b = 0.95

  density_r_2 = 15
  albedo_r_2 = 0.95

  density_g_2 = 100
  albedo_g_2 = 0.001

  density_b_2 = 30
  albedo_b_2 = 0.001

  opts = {
    'scene_id': 'smoke_two_ss', 
    'g' : -0.4, 
    'l' : 0.5,
    'gamma': 0.15,
    's_r' : density_r * albedo_r, 
    's_g' : density_g * albedo_g, 
    's_b' : density_b * albedo_b, 
    'a_r' : density_r * (1.0 - albedo_r),
    'a_g' : density_g * (1.0 - albedo_g),
    'a_b' : density_b * (1.0 - albedo_b),
    's_r_2' : density_r_2 * albedo_r_2, 
    's_g_2' : density_g_2 * albedo_g_2, 
    's_b_2' : density_b_2 * albedo_b_2, 
    'a_r_2' : density_r_2 * (1.0 - albedo_r_2), 
    'a_g_2' : density_g_2 * (1.0 - albedo_g_2),
    'a_b_2' : density_b_2 * (1.0 - albedo_b_2),
    'noise_x': noise_x,
    'noise_y': noise_y,
    'noise_z': noise_z,
  }

  args.custom_formatting.append(opts)
  args.threads = 8
  return args

def setup_colored_smoke_three(args):
  args.xres = 800
  args.yres = 1600
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure5/smoke_template.pbrt"]
  args.integrators = ["nullpath"]
  args.integrator_flags = {
      "spectral" : [
           #{"flags": [COMMON_FLAGS.enable_split_nee,
           #           COMMON_FLAGS.single_scatter,
           #           COMMON_FLAGS.ignore_no_scatter],
           #  "samples": 3000},
           #{"flags": [COMMON_FLAGS.enable_split_nee,
           #          COMMON_FLAGS.single_scatter,
           #          COMMON_FLAGS.ignore_no_scatter],
           #  "samples": 28},
            #{"flags": [COMMON_FLAGS.enable_naive_delta,
            #           COMMON_FLAGS.enable_split_nee,
            #           COMMON_FLAGS.single_scatter,
            #           COMMON_FLAGS.ignore_no_scatter],
            #"samples": 31},
       ],
      "nullpath" : [
        #{"flags" : [COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.single_scatter,
        #            COMMON_FLAGS.ignore_no_scatter],
        #    "samples": 512},
        #{"flags" : [COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.single_scatter,
        #            COMMON_FLAGS.ignore_no_scatter],
        #    "samples" : 28},
        {"flags" : [COMMON_FLAGS.enable_delta,
                   COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.single_scatter,
                    COMMON_FLAGS.ignore_no_scatter],
            "samples" : 30},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  noise_x = 0.25
  noise_y = 3
  noise_z = 0.25

  density_r = 30
  albedo_r = 0.001

  density_g = 35
  albedo_g = 0.95

  density_b = 15
  albedo_b = 0.95

  density_r_2 = 15
  albedo_r_2 = 0.95

  density_g_2 = 15
  albedo_g_2 = 0.001

  density_b_2 = 30
  albedo_b_2 = 0.001

  opts = {
    'scene_id': 'smoke_three_ss', 
    'g' : -0.4, 
    'l' : 0.5,
    'gamma': 0.15,
    's_r' : density_r * albedo_r, 
    's_g' : density_g * albedo_g, 
    's_b' : density_b * albedo_b, 
    'a_r' : density_r * (1.0 - albedo_r),
    'a_g' : density_g * (1.0 - albedo_g),
    'a_b' : density_b * (1.0 - albedo_b),
    's_r_2' : density_r_2 * albedo_r_2, 
    's_g_2' : density_g_2 * albedo_g_2, 
    's_b_2' : density_b_2 * albedo_b_2, 
    'a_r_2' : density_r_2 * (1.0 - albedo_r_2), 
    'a_g_2' : density_g_2 * (1.0 - albedo_g_2),
    'a_b_2' : density_b_2 * (1.0 - albedo_b_2),
    'noise_x': noise_x,
    'noise_y': noise_y,
    'noise_z': noise_z,
  }

  args.custom_formatting.append(opts)
  args.threads = 8
  return args

def setup_colored_smoke_one_ms(args):
  args.xres = 800
  args.yres = 1600
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure5/smoke_template.pbrt"]
  args.integrators = ["spectral"]
  args.integrator_flags = {
      "spectral" : [
        {"flags": [COMMON_FLAGS.enable_split_nee,
                     COMMON_FLAGS.ignore_no_scatter],
         "samples": 3000},
        #{"flags": [COMMON_FLAGS.enable_split_nee,
        #             COMMON_FLAGS.enable_naive_delta,
        #             COMMON_FLAGS.ignore_no_scatter],
        # "samples": 32},
        #{"flags": [COMMON_FLAGS.enable_split_nee,
        #         COMMON_FLAGS.ignore_no_scatter],
        # "samples": 16}
      ],
      "nullpath" : [
        {"flags" : [COMMON_FLAGS.enable_delta,
                    COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.ignore_no_scatter],
         "samples": 32},
      #  {"flags" : [COMMON_FLAGS.enable_split_nee,
      #              COMMON_FLAGS.ignore_no_scatter],
      #   "samples" : 31},
      #{"flags" : [COMMON_FLAGS.enable_split_nee,
      #            COMMON_FLAGS.ignore_no_scatter],
      # "samples" : 512},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  noise_x = 0.5
  noise_y = 4
  noise_z = 0.5

  density_r = 30
  albedo_r = 0.001

  density_g = 100
  albedo_g = 0.001

  density_b = 15
  albedo_b = 0.95

  density_r_2 = 15
  albedo_r_2 = 0.95

  density_g_2 = 100
  albedo_g_2 = 0.001

  density_b_2 = 30
  albedo_b_2 = 0.001

  opts = {
    'scene_id': 'smoke_one_ms', 
    'g' : -0.4, 
    'l' : 0.75,
    'gamma': 0.15,
    's_r' : density_r * albedo_r, 
    's_g' : density_g * albedo_g, 
    's_b' : density_b * albedo_b, 
    'a_r' : density_r * (1.0 - albedo_r),
    'a_g' : density_g * (1.0 - albedo_g),
    'a_b' : density_b * (1.0 - albedo_b),
    's_r_2' : density_r_2 * albedo_r_2, 
    's_g_2' : density_g_2 * albedo_g_2, 
    's_b_2' : density_b_2 * albedo_b_2, 
    'a_r_2' : density_r_2 * (1.0 - albedo_r_2), 
    'a_g_2' : density_g_2 * (1.0 - albedo_g_2),
    'a_b_2' : density_b_2 * (1.0 - albedo_b_2),
    'noise_x': noise_x,
    'noise_y': noise_y,
    'noise_z': noise_z,
  }

  args.custom_formatting.append(opts)
  args.threads = 16
  return args

def setup_colored_smoke_one_new(args):
  args.xres = 800
  args.yres = 1600
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure5/smoke_template.pbrt"]
  args.integrators = ["nullpath"] #, "nullpath"]
  args.integrator_flags = {
      "spectral" : [
          #{"flags": [COMMON_FLAGS.enable_split_nee,
          #           COMMON_FLAGS.single_scatter,
          #           COMMON_FLAGS.ignore_no_scatter],
          # "samples": 3000},
          #{"flags": [COMMON_FLAGS.enable_split_nee,
          #           COMMON_FLAGS.enable_naive_delta,
          #           COMMON_FLAGS.single_scatter,
          #           COMMON_FLAGS.ignore_no_scatter],
          #           "samples": 32},
         #{"flags": [COMMON_FLAGS.enable_split_nee,
         #            COMMON_FLAGS.single_scatter,
         #            COMMON_FLAGS.ignore_no_scatter],
         #            "samples": 16},
       ],
      "nullpath" : [
        {"flags" : [COMMON_FLAGS.enable_delta,
                    COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.single_scatter,
                    COMMON_FLAGS.ignore_no_scatter],
            "samples" : 30},
        #{"flags" : [COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.single_scatter,
        #            COMMON_FLAGS.ignore_no_scatter],
        #    "samples" : 28},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  noise_x = 0.5
  noise_y = 4
  noise_z = 0.5

  density_r = 30
  albedo_r = 0.95

  density_g = 100
  albedo_g = 0.001

  density_b = 15
  albedo_b = 0.001

  density_r_2 = 15
  albedo_r_2 = 0.001

  density_g_2 = 100
  albedo_g_2 = 0.001

  density_b_2 = 30
  albedo_b_2 = 0.95

  opts = {
    'scene_id': 'smoke_one_ms', 
    'g' : -0.4, 
    'l' : 0.75,
    'gamma': 0.15,
    's_r' : density_r * albedo_r, 
    's_g' : density_g * albedo_g, 
    's_b' : density_b * albedo_b, 
    'a_r' : density_r * (1.0 - albedo_r),
    'a_g' : density_g * (1.0 - albedo_g),
    'a_b' : density_b * (1.0 - albedo_b),
    's_r_2' : density_r_2 * albedo_r_2, 
    's_g_2' : density_g_2 * albedo_g_2, 
    's_b_2' : density_b_2 * albedo_b_2, 
    'a_r_2' : density_r_2 * (1.0 - albedo_r_2), 
    'a_g_2' : density_g_2 * (1.0 - albedo_g_2),
    'a_b_2' : density_b_2 * (1.0 - albedo_b_2),
    'noise_x': noise_x,
    'noise_y': noise_y,
    'noise_z': noise_z,
  }

  args.custom_formatting.append(opts)
  args.threads = 8
  return args

def setup_colored_smoke_two_ms(args):
  args.xres = 800
  args.yres = 1600
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure5/smoke_template.pbrt"]
  args.integrators = ["nullpath"]
  args.integrator_flags = {
      "spectral" : [
         #{"flags": [COMMON_FLAGS.enable_split_nee,
         #            COMMON_FLAGS.ignore_no_scatter], 
         #  "samples": 3000},
          {"flags": [COMMON_FLAGS.enable_split_nee,
                     COMMON_FLAGS.enable_naive_delta,
                     COMMON_FLAGS.ignore_no_scatter], 
           "samples": 32},
          #{"flags": [COMMON_FLAGS.enable_split_nee,
          #           COMMON_FLAGS.ignore_no_scatter],
          #  "samples": 17},
       ],
      "nullpath" : [
        #{"flags" : [COMMON_FLAGS.enable_delta,
        #            COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.ignore_no_scatter],
        # "samples": 34},
        {"flags" : [COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.ignore_no_scatter],
         "samples": 29},
        #{"flags" : [COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.ignore_no_scatter],
        #    "samples" : 512},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  noise_x = 0.5
  noise_y = 2
  noise_z = 0.5

  density_r = 30
  albedo_r = 0.001

  density_g = 100
  albedo_g = 0.95

  density_b = 15
  albedo_b = 0.95

  density_r_2 = 15
  albedo_r_2 = 0.95

  density_g_2 = 100
  albedo_g_2 = 0.001

  density_b_2 = 30
  albedo_b_2 = 0.001

  opts = {
    'scene_id': 'smoke_two_ms', 
    'g' : -0.4, 
    'l' : 0.75,
    'gamma': 0.15,
    's_r' : density_r * albedo_r, 
    's_g' : density_g * albedo_g, 
    's_b' : density_b * albedo_b, 
    'a_r' : density_r * (1.0 - albedo_r),
    'a_g' : density_g * (1.0 - albedo_g),
    'a_b' : density_b * (1.0 - albedo_b),
    's_r_2' : density_r_2 * albedo_r_2, 
    's_g_2' : density_g_2 * albedo_g_2, 
    's_b_2' : density_b_2 * albedo_b_2, 
    'a_r_2' : density_r_2 * (1.0 - albedo_r_2), 
    'a_g_2' : density_g_2 * (1.0 - albedo_g_2),
    'a_b_2' : density_b_2 * (1.0 - albedo_b_2),
    'noise_x': noise_x,
    'noise_y': noise_y,
    'noise_z': noise_z,
  }

  args.custom_formatting.append(opts)
  args.threads = 8
  return args

def setup_colored_smoke_three_ms(args):
  args.xres = 800
  args.yres = 1600
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure5/smoke_template.pbrt"]
  args.integrators = ["nullpath"]
  args.integrator_flags = {
      "spectral" : [
        #{"flags": [COMMON_FLAGS.enable_split_nee,
        #           COMMON_FLAGS.ignore_no_scatter],
        #  "samples": 3000}],
        {"flags": [COMMON_FLAGS.enable_split_nee,
                   COMMON_FLAGS.ignore_no_scatter],
         "samples": 29},
        #{"flags": [COMMON_FLAGS.enable_naive_delta,
        #           COMMON_FLAGS.enable_split_nee,
        #           COMMON_FLAGS.ignore_no_scatter],
        # "samples": 32}
        ],
      "nullpath" : [
        #{"flags" : [COMMON_FLAGS.enable_delta,
        #            COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.ignore_no_scatter],
        # "samples": 31},
        {"flags" : [COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.ignore_no_scatter],
         "samples": 27},
        #{"flags" : [COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.ignore_no_scatter],
        # "samples" : 512},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  noise_x = 0.25
  noise_y = 3
  noise_z = 0.25

  density_r = 30
  albedo_r = 0.001

  density_g = 35
  albedo_g = 0.95

  density_b = 15
  albedo_b = 0.95

  density_r_2 = 15
  albedo_r_2 = 0.95

  density_g_2 = 15
  albedo_g_2 = 0.001

  density_b_2 = 30
  albedo_b_2 = 0.001

  opts = {
    'scene_id': 'smoke_three_ms', 
    'g' : -0.4, 
    'l' : 0.5,
    'gamma': 0.15,
    's_r' : density_r * albedo_r, 
    's_g' : density_g * albedo_g, 
    's_b' : density_b * albedo_b, 
    'a_r' : density_r * (1.0 - albedo_r),
    'a_g' : density_g * (1.0 - albedo_g),
    'a_b' : density_b * (1.0 - albedo_b),
    's_r_2' : density_r_2 * albedo_r_2, 
    's_g_2' : density_g_2 * albedo_g_2, 
    's_b_2' : density_b_2 * albedo_b_2, 
    'a_r_2' : density_r_2 * (1.0 - albedo_r_2), 
    'a_g_2' : density_g_2 * (1.0 - albedo_g_2),
    'a_b_2' : density_b_2 * (1.0 - albedo_b_2),
    'noise_x': noise_x,
    'noise_y': noise_y,
    'noise_z': noise_z,
  }

  args.custom_formatting.append(opts)
  args.threads = 16
  return args


def setup_colored_smoke_teaser(args):
  args.xres = 2400
  args.yres = 1600
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/teaser/smoke_template.pbrt"]
  args.integrators = ["spectral"] #, "spectral"]
  args.integrator_flags = {
      "spectral" : [{
          "flags": [COMMON_FLAGS.enable_naive_delta,
                    COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.single_scatter,
                    COMMON_FLAGS.ignore_no_scatter],
          "samples": 3000},
        #{"flags": [COMMON_FLAGS.enable_split_nee,
         #          COMMON_FLAGS.single_scatter,
         #          COMMON_FLAGS.ignore_no_scatter],
         # "samples": 8},
         #{"flags": [COMMON_FLAGS.enable_naive_delta,
         #           COMMON_FLAGS.enable_split_nee,
         #           COMMON_FLAGS.single_scatter,
         #           COMMON_FLAGS.ignore_no_scatter],
         # "samples": 16}
      ],
      "nullpath" : [
        {"flags" : [COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.single_scatter,
                    COMMON_FLAGS.ignore_no_scatter],
         "samples" : 15},
        {"flags" : [COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.single_scatter,
                    COMMON_FLAGS.ignore_no_scatter],
         "samples" : 128}],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  noise_x = 0.5
  noise_y = 2
  noise_z = 0.5

  density_r = 15
  albedo_r = 0.001

  density_g = 100
  albedo_g = 0.001

  density_b = 30
  albedo_b = 0.95

  density_r_2 = 30
  albedo_r_2 = 0.95

  density_g_2 = 100
  albedo_g_2 = 0.001

  density_b_2 = 15
  albedo_b_2 = 0.001

  opts = {
    'scene_id': 'teaser', 
    'g' : -0.4, 
    'l' : 0.5,
    'gamma': 0.15,
    's_r' : density_r * albedo_r, 
    's_g' : density_g * albedo_g, 
    's_b' : density_b * albedo_b, 
    'a_r' : density_r * (1.0 - albedo_r),
    'a_g' : density_g * (1.0 - albedo_g),
    'a_b' : density_b * (1.0 - albedo_b),
    's_r_2' : density_r_2 * albedo_r_2, 
    's_g_2' : density_g_2 * albedo_g_2, 
    's_b_2' : density_b_2 * albedo_b_2, 
    'a_r_2' : density_r_2 * (1.0 - albedo_r_2), 
    'a_g_2' : density_g_2 * (1.0 - albedo_g_2),
    'a_b_2' : density_b_2 * (1.0 - albedo_b_2),
    'noise_x': noise_x,
    'noise_y': noise_y,
    'noise_z': noise_z,
  }

  args.custom_formatting.append(opts)
  args.threads = 16
  return args

def setup_ea_cloud(args):
  args.xres = 1024
  args.yres = 768
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure6/cloud_template.pbrt"]
  args.integrators = ["nullpath"]
  args.integrator_flags = {
    #"spectral": [
      #{"flags" : [COMMON_FLAGS.enable_split_nee,
      #            COMMON_FLAGS.ignore_no_scatter,
      #            COMMON_FLAGS.single_scatter],
      # "samples" : 32}],
      "nullpath": [
          {"flags" : [COMMON_FLAGS.enable_only_nee,
                      COMMON_FLAGS.ignore_no_scatter,
                      COMMON_FLAGS.single_scatter],
            "samples" : 30},
          #{"flags" : [COMMON_FLAGS.enable_only_ea,
          #            COMMON_FLAGS.ignore_no_scatter,
          #            COMMON_FLAGS.single_scatter],
          #  "samples" : 27},
          #{"flags" : [COMMON_FLAGS.enable_dist_mis,
          #            COMMON_FLAGS.ignore_no_scatter,
          #            COMMON_FLAGS.single_scatter],
          #  "samples" : 20
          #},
          #{"flags": [COMMON_FLAGS.enable_split_nee,
          #           COMMON_FLAGS.ignore_no_scatter,
          #           COMMON_FLAGS.single_scatter],
          #   "samples": 25
          #},
          #{"flags": [COMMON_FLAGS.enable_dist_mis,
          #           COMMON_FLAGS.ignore_no_scatter,
          #           COMMON_FLAGS.single_scatter],
          #   "samples": 3000
          #}
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  g = 0
  albedo = 1
  d = 20
  s = d * float(albedo)
  a = d - s

  r1 = 0.1
  translate1 = '0 2 0'
  l1 = '150 150 150'

  r2 = 0.0001
  translate2 = '-0.2 -0.7 0.2'
  l2 = '5000000 0 0'

  r3 = 0.0001
  translate3 = '-0.3 -0.2 0.35'
  l3 = '0 0 5000000'

  opts = {
    'scene_id' : 'ea_cloud', 
    'g' : g, 
    'r1' : r1,
    'r2' : r2,
    'r3' : r3,
    's' : s,
    'a' : a,
    'l1' : l1,
    'l2' : l2,
    'l3' : l3,
    'translate1': translate1,
    'translate2': translate2,
    'translate3': translate3,
    'gamma': 1,
    'useLinearRamp': 'false',
    'rampScale': -0.04,
    'rampDimension': 0,
    'rampOrigin': 0.4,
  }
  args.custom_formatting.append(opts)
  args.threads = 8
  return args

def setup_dense_dragon(args):
  args.xres = 1024
  args.yres = 768
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/figure2/scene_template.pbrt"]
  args.integrators = ["nullpath"]
  args.integrator_flags = {
    "spectral": [
        #{ "flags" : [], "samples" : 30000},
        #{ "flags" : [], "samples" : 36 },
        #{ "flags" : [], "samples": 72 },
        #{ "flags" : [], "samples": 144 },
        #{ "flags" : [], "samples": 288 },
        #{ "flags" : [], "samples": 576 },
        #{ "flags" : [COMMON_FLAGS.enable_naive_delta], "samples" : 90 }
        #{ "flags" : [COMMON_FLAGS.enable_naive_delta], "samples" : 180 }
        #{ "flags" : [COMMON_FLAGS.enable_naive_delta], "samples" : 360 }
        #{ "flags" : [COMMON_FLAGS.enable_naive_delta], "samples" : 720 }
        #{ "flags" : [COMMON_FLAGS.enable_naive_delta], "samples" : 1440 }
    ],
    "nullpath": [
        { "flags": [], "samples": 64 },
        #{ "flags": [], "samples": 128 },
        #{ "flags": [], "samples": 256 },
        #{ "flags": [], "samples": 512 },
        #{ "flags": [], "samples": 1024 },
    ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  scale = 10
  # Medium Params 1
  r1_density = 60 * scale
  g1_density = 5 * scale
  b1_density = 6 * scale
  
  r1_albedo = 0.001
  g1_albedo = 0.0075
  b1_albedo = 0.001

  r1_null = 1 * scale
  g1_null = 1 * scale
  b1_null = 1 * scale

  # Medium Params 2
  r2_density = 2 * scale
  g2_density = 2.5 * scale
  b2_density = 3 * scale

  r2_albedo = 0.1
  g2_albedo = 0.7
  b2_albedo = 0.35

  r2_null = 1 * scale
  g2_null = 1 * scale
  b2_null = 1 * scale

  opts = {
    'scene_id' : 'dense_dragon',
    'eta': 1.2,
    'g' : 0,
    'w' : 0.5,
    'num_octaves': 8,
    'noise': 3,
    'noise2': 10,
    'noise3': 37,
    'noise4': 0.8,
    's1_r' : r1_density * r1_albedo,
    'a1_r' : r1_density * (1.0 - r1_albedo),
    'n1_r' : r1_null,
    's1_g' : g1_density * g1_albedo,
    'a1_g' : g1_density * (1.0 - g1_albedo),
    'n1_g' : g1_null,
    's1_b' : b1_density * b1_albedo,
    'a1_b' : b1_density * (1.0 - b1_albedo),
    'n1_b' : b1_null,
    's2_r' : r2_density * r2_albedo,
    'a2_r' : r2_density * (1.0 - r2_albedo),
    'n2_r' : r2_null,
    's2_g' : g2_density * g2_albedo,
    'a2_g' : g2_density * (1.0 - g2_albedo),
    'n2_g' : g2_null,
    's2_b' : b2_density * b2_albedo,
    'a2_b' : b2_density * (1.0 - b2_albedo),
    'n2_b' : b2_null,
  }
  args.custom_formatting.append(opts)
  args.threads = 16
  return args

def setup_breakfast(args):
  args.xres = 1024
  args.yres = 768
  args.scenes = ["/home/johnn/Dropbox/research/volumemis/figures/figure7/breakfast.pbrt"]
  args.integrators = ["spectral", "nullpath"]
  args.integrator_flags = {
    "spectral": [
        {
            "flags": [
                 COMMON_FLAGS.enable_split_nee,
                 #COMMON_FLAGS.enable_surface_nee,
            ],
            "samples": 128
        }
    ],
    "nullpath": [
        {
            "flags": [
                COMMON_FLAGS.enable_split_nee,
                #COMMON_FLAGS.enable_surface_nee,
            ],
            "samples": 114
        },
        #{
        #    "flags": [
        #        #COMMON_FLAGS.enable_split_nee,
        #        #COMMON_FLAGS.enable_delta,
        #        #COMMON_FLAGS.enable_surface_nee,
        #     ],
        #    "samples": 10
        #}
    ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  # Medium Params 1
  density = 0.35
  albedo = 0.3
  null = 0.2

  opts = {
    'scene_id' : 'breakfast',
    'g' : 0, 
    'gamma': 1,
    's' : density * albedo,
    'a' : density * (1.0 - albedo),
    'n' : null,
  }
  args.custom_formatting.append(opts)
  args.threads = 16
  return args

def setup_cornell(args):
  args.xres = 1024
  args.yres = 1024
  args.scenes = ["/home/bailey/Dropbox/research/volumemis/figures/cornell/cornell.pbrt"]
  args.integrators = ["nullpath"]
  args.integrator_flags = {
    "nullpath": [
        { 
            "flags": [
                COMMON_FLAGS.maxDepth(4),
                COMMON_FLAGS.enable_split_nee,
                COMMON_FLAGS.enable_surface_nee,
                COMMON_FLAGS.ignore_no_scatter,
            ], 
            "samples": 512,
        },
     ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []

  opts = {
    'scene_id' : 'dense_dragon',
  }
  args.custom_formatting.append(opts)
  args.threads = 16
  return args

def setup_furnace_path_view(args):
  args.xres = 20
  args.yres = 20
  args.scenes = ["../scenes/final/furnace/furnace_with_sphere_and_env_light.pbrt"]
  args.integrators = ["nullpath"]
  numEvents = 1
  numSamples = 1
  args.integrator_flags = {
      "nullpath" : [
        {"flags": [COMMON_FLAGS.enable_nee,
                   COMMON_FLAGS.print_path], "samples": numSamples},
      ],
  }
  
  args.custom_scene_file = True
  args.custom_formatting = []
  s = 1
  a = 0
  g = 0
  args.custom_formatting.append({
      'scene_id' : 'albedo={0}_density={1}='.format(s / float(s + a), s +  a), 
      'g' : g, 
      's' : s,
      'a' : a,
  })
  args.threads = 1
  return args

def setup_furnace_env_light(args):
  args.xres = 1
  args.yres = 1
  #args.scenes = ["../scenes/final/furnace/furnace_with_env_light.pbrt"]
  args.scenes = ["../scenes/final/furnace/furnace_with_sphere_and_env_light.pbrt"]
  args.integrators = ["nullpath"]
  numSamples = 1000000
  args.integrator_flags = {
      "nullpath" : [
        {"flags": [], "samples": numSamples},
        {"flags": [COMMON_FLAGS.enable_nee],
          "samples": numSamples},
        {"flags": [COMMON_FLAGS.enable_only_nee],
          "samples": numSamples},
        {"flags": [COMMON_FLAGS.enable_split_nee],
          "samples": numSamples},
      ],
  }
  
  args.custom_scene_file = True
  args.custom_formatting = []
  s = 20
  a = 0
  g = 0
  args.custom_formatting.append({
    'scene_id' : 'albedo={0}_density={1}_xres={2}'.format(s / float(s + a), s +  a, args.xres), 
    'g' : g, 
    's' : s,
    'a' : a,
  })
  args.threads = 16
  return args

def setup_ea_furnace(args):
  args.xres = 100
  args.yres = 100
  args.scenes = [
   "../scenes/final/infinite/infinite_ea_template.pbrt"
  ]
  args.integrators = ["nullpath"]
  numSamples = 100
  args.integrator_flags = {
      "nullpath" : [
        {"flags": [COMMON_FLAGS.enable_nee,
                   COMMON_FLAGS.viewNumScatterEvents(2)], 
          "samples": numSamples},
        {"flags": [COMMON_FLAGS.enable_only_nee,
                   COMMON_FLAGS.viewNumScatterEvents(2)], 
         "samples": numSamples},
        {"flags": [COMMON_FLAGS.enable_only_ea,
                   COMMON_FLAGS.viewNumScatterEvents(2)], 
         "samples": numSamples},
        {"flags": [COMMON_FLAGS.enable_dist_mis,
                    COMMON_FLAGS.viewNumScatterEvents(2)], 
         "samples": numSamples},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []
  s = 1
  a = 0
  n = 10
  g = 0
  L = 100000
  args.custom_formatting.append({
      'scene_id' : 'albedo={0}_density={1}_n={2}_xres={3}'.format(s / float(s + a), s +  a, n, args.xres), 
      'g' : g, 
      's' : s,
      'a' : a,
      'n' : n,
      'r' : 0.05,
      'translate': '-5 0 0',
      'L' : L,
  })
  args.threads = 1
  return args

def setup_furnace(args):
  args.xres = 1000
  args.yres = 1000
  #args.scenes = [
  # "../scenes/final/furnace/furnace_with_no_emission_light.pbrt"
  #]
  args.scenes = [
   "../scenes/final/furnace/furnace_with_two_no_emission_lights.pbrt"
  ]
  #args.scenes = ["../scenes/final/furnace/furnace_template.pbrt"]
  args.integrators = ["nullpath"]
  numEvents = 0
  numSamples = 1
  args.integrator_flags = {
      "nullpath" : [
        {"flags": [COMMON_FLAGS.enable_nee,
                   COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(2)],
                   "samples": numSamples},
        {"flags": [COMMON_FLAGS.enable_nee,
                   COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(2)],
                   "samples": numSamples},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(1)],
          "samples": numSamples},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(2)],
          "samples": numSamples},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(3)],
          "samples": numSamples},
        {"flags": [COMMON_FLAGS.enable_nee,
                   COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(numEvents)], "samples": numSamples},
        {"flags": [COMMON_FLAGS.enable_only_nee, 
                   COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.ignore_no_scatter], "samples": numSamples},
        {"flags": [COMMON_FLAGS.enable_split_nee, 
                   COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.ignore_no_scatter], "samples": numSamples},
      ],
  }
  
  args.custom_scene_file = True
  args.custom_formatting = []
  s = 1
  a = 1
  g = 0.8
  args.custom_formatting.append({
      'scene_id' : 'albedo={0}_density={1}_xres={2}'.format(s / float(s + a), s +  a, args.xres), 
      'g' : g, 
      's' : s,
      'a' : a,
  })
  args.threads = 3
  return args

def setup_infinite(args):
  args.scenes = ["../scenes/final/infinite/infinite_template.pbrt"]
  args.integrators = ["nullpath"]
  args.integrator_flags = {
      "nullpath" : [
        {"flags": [COMMON_FLAGS.constantEmission(1)],
            "samples": 1},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []
  s = 20
  a = 20
  g = 0.9
  args.custom_formatting.append({
      'scene_id' : 'g={0}_s={1}'.format(g,s), 
      'g' : g, 
      's' : s,
      'a' : a,
  })
  args.threads = 1
  return args

def setup_infinite_avg(args):
  args.xres = 1
  args.yres = 1
  args.scenes = ["../scenes/final/infinite/infinite_template.pbrt"]
  args.integrators = ["nullpath"]
  args.integrator_flags = {
      "nullpath" : [
        {"flags": [COMMON_FLAGS.constantEmission(1)],
            "samples": 480000},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(0)],
            "samples": 480000},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(1)],
            "samples": 480000},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(2)],
            "samples": 480000},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(3)],
            "samples": 480000},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(4)],
            "samples": 480000},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.enable_nee],
            "samples": 480000},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(0),
                   COMMON_FLAGS.enable_nee],
            "samples": 480000},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(1),
                   COMMON_FLAGS.enable_nee],
            "samples": 480000},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(2),
                   COMMON_FLAGS.enable_nee],
            "samples": 480000},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(3),
                   COMMON_FLAGS.enable_nee],
            "samples": 480000},
        {"flags": [COMMON_FLAGS.constantEmission(1),
                   COMMON_FLAGS.viewNumScatterEvents(4),
                   COMMON_FLAGS.enable_nee],
            "samples": 480000},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []
  g = 0.9
  s = 20
  a = 20
  args.custom_formatting.append({
      'scene_id' : 'avg_g={0}_s={1}'.format(g, s), 
      'g' : g, 
      's' : s,
      'a' : a,
  })
  args.threads = 1
  return args

def setup_two_light_comparison(args):
  args.xres = 500
  args.yres = 500
  args.scenes = ["../scenes/final/veach_bunny/two_light_bunny_template.pbrt"]
  args.integrators = ["nullpath"]
  args.integrator_flags = {
      "nullpath" : [
        #{"flags" : [COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.show_nee_mis_weights,
              #      COMMON_FLAGS.single_scatter,
        #            COMMON_FLAGS.ignore_no_scatter],
        #"samples" : 16},
        #{"flags" : [COMMON_FLAGS.enable_split_nee,
        #            COMMON_FLAGS.split_nee_mis_values,
        #            COMMON_FLAGS.ignore_no_scatter,
        #             COMMON_FLAGS.single_scatter],
        #"samples" : 1},
        {"flags" : [
                     COMMON_FLAGS.enable_naive,
                     COMMON_FLAGS.ignore_no_scatter],
        "samples" : 2048},
        {"flags" : [COMMON_FLAGS.enable_nee,
                    COMMON_FLAGS.ignore_no_scatter],
        "samples" : 2048},
        {"flags" : [COMMON_FLAGS.enable_only_nee,
                    COMMON_FLAGS.ignore_no_scatter],
        "samples" : 2048},
        {"flags" : [COMMON_FLAGS.enable_split_nee,
                    COMMON_FLAGS.ignore_no_scatter],
        "samples" : 2048},
      ],
  }
  args.custom_scene_file = True
  args.custom_formatting = []
 
  
  albedo = 0.75
  L1 = 15
  L2 = 5
  g_opts = [0, 0.9, 0.99]
  density_opts = [1, 10, 20]  
  r1 = 0.25
  r2 = 1

  for g, d in itertools.product(g_opts, density_opts):
    s = d * float(albedo)
    a = d - s
    translate1 = '0 {0} 0'.format(-r1 - 1.1)
    translate2 = '{0} 0 0'.format(-r2 - 5) 
    l1 = L1 /float(pow(r1, 2.0))
    l2 = L2 /float(pow(r2, 2.0))
    opts = {
        'scene_id' : 'g={0}_d={1}'.format(g,d), 
        'g' : g, 
        'r1' : r1,
        'r2' : r2,
        's' : s,
        'a' : a,
        'l1' : l1,
        'l2' : l2,
        'translate1': translate1,
        'translate2': translate2,
    }
    args.custom_formatting.append(opts)
  args.threads = 2
  return args

# Create scene files generates temporary, custom scene files:
#   - Generates a template scene file
#   - Creates "j" temporary scene files (one for each job)
def create_custom_scene_files(args):
  pbrt_files_to_run = []
  for scene, scene_id in create_custom_scene_template(args.scene, args.custom_formatting):
      # For each integrator, create a scene file for each job, and each set 
      # of flags.
      for integrator in args.integrators:
        
        # If there are no flags, create an empty flag list for the integrator.
        int_flags = args.integrator_flags
        flag_lists = [[]] if integrator not in int_flags else int_flags[integrator]
        
        for f in flag_lists:
            if type(f) is dict:
              flag_list = f["flags"] 
              samples = args.samples if "samples" not in f else f["samples"]
            else:
              flag_list = f
              samples = args.samples

            # Get scene file info.
            input_file, scene_ext  = os.path.splitext(args.scene)
            scene_dir = os.path.abspath(os.path.dirname(input_file))
            scene_file = os.path.basename(input_file)
            
            # Determine output file name.
            flags = ""
            identifier = "integrator={0}".format(integrator)
            for flag in flag_list:
              if "name" not in flag:
                flag["name"] = None

              #Add the flag to the file name.
              flagDisplay = flag["name"] if "display" not in flag else flag["display"]
              identifier += "_" + str(flagDisplay) + "_" 
             
              if flag["name"] == None:
                continue

              # Format the flag as a PBRT Integrator flag.
              flags += "\n\t\"" + flag["type"] + " "+ str(flag["name"]) + "\" " + str(flag["value"])
           
            for j in range(args.jobs):
              # if there is only one job do not include the job #
              if args.jobs == 1:
                j = None

              # Determine temporary scene file name.
              input_file = os.path.join(scene_dir, 
                get_scene_filename(scene_file, scene_id, identifier, samples, j, scene_ext))
              
              pbrt_files_to_run.append(input_file)

              # Format a temporary scene file
              formatted_scene = scene.format(
                args.xres, args.yres,
                os.path.join(args.out_dir,
                             get_image_filename(scene_file, scene_id, identifier, 
                                                samples, j)),
                args.scale,
                samples,
                integrator,
                flags)

              # Save
              f = open(input_file, "w")
              f.write(formatted_scene)
              f.close()
  
  return pbrt_files_to_run

def create_scene_files(args):
  scene = create_scene_template(args.scene)
  pbrt_files_to_run = []

  # For each integrator, create a scene file for each job, and each set 
  # of flags.
  for integrator in args.integrators:
    
    # If there are no flags, create an empty flag list for the integrator.
    int_flags = args.integrator_flags
    flag_lists = [[]] if integrator not in int_flags else int_flags[integrator]
    
    for f in flag_lists:
      if type(f) is dict:
        flag_list = f["flags"] 
        samples = args.samples if "samples" not in f else f["samples"]
      else:
        flag_list = f
        samples = args.samples
        
      # Get scene file info.
      input_file, scene_ext  = os.path.splitext(args.scene)
      scene_dir = os.path.abspath(os.path.dirname(input_file))
      scene_file = os.path.basename(input_file)
        
      # Determine output file name.
      flags = ""
      identifier = integrator
      for flag in flag_list:
        identifier += "_" + str(flag["name"]) + "_" 
        flags += "\n\t\"" + flag["type"] + " "+ str(flag["name"]) + "\" " + str(flag["value"])
       
      for j in range(args.jobs):
        # Determine temporary scene file name.
        input_file = os.path.join(scene_dir,
                get_scene_filename(scene_file, identifier, samples, j, scene_ext))
          
        pbrt_files_to_run.append(input_file)
        # Format a temporary scene file
        formatted_scene = scene.format(
            args.xres, args.yres,
            os.path.join(args.out_dir,
                     get_image_filename(scene_file, identifier, 
                                        samples, j)),
            samples,
            integrator,
            flags)
        
        # Save
        f = open(input_file, "w")
        f.write(formatted_scene)
        f.close()
  
  return pbrt_files_to_run

# Create a template scene file that can be adjusted with
# formatting. The params added are the following:
#
# 0.) x-resolution
# 1.) y-resolution
# 2.) output filename
# 3.) samples per pixel
# 4.) integrator
# 5.) integrator flags
#
def create_scene_template(scene_file): 
  custom_scene_data = "Film \"image\"\
\n\t\"integer xresolution\" [{0}]\
\n\t\"integer yresolution\" [{1}]\
\n\t\"string filename\" \"{2}\"\n\
Sampler \"random\" \"integer pixelsamples\" [{3}]\n\
Integrator \"{4}\" {5}\
"
  with open(scene_file, 'r') as input_scene:
    scene_data = input_scene.read()
  return custom_scene_data + "\n" + scene_data

# Create a custom template scene file that can be adjusted with
# formatting. The params added are the following:
#
# 0.) x-resolution
# 1.) y-resolution
# 2.) output filename
# 3.) samples per pixel
# 4.) integrator
# 5.) integrator flags
#
def create_custom_scene_template(scene_file, custom_formats): 
  custom_film = "Film \"image\"\
\n\t\"integer xresolution\" [{0}]\
\n\t\"integer yresolution\" [{1}]\
\n\t\"string filename\" \"{2}\"\
\n\t\"float scale\" {3}"
  custom_sampler_integrator = "\
Sampler \"random\" \"integer pixelsamples\" [{4}]\n\
Integrator \"{5}\" {6}\
"
  with open(scene_file, 'r') as input_scene:
    scene_data = input_scene.read()
  scenes = []
  for custom_format in custom_formats:
    if 'cropwindow' not in custom_format:
        custom_format['cropwindow'] = '0 1.0 0 1.0'
    crop_window = '\n\t\"float cropwindow\" [' + custom_format['cropwindow'] + ']\n'
    custom_scene_data = custom_film + crop_window + custom_sampler_integrator
    scenes.append((custom_scene_data + "\n" + scene_data.format(**custom_format), custom_format['scene_id']))
  return scenes

def get_image_filename(scene_file, scene_id, identifier, samples, job, image_ext = ".exr"):
  if job:
    return "scene={0}_{1}_{2}samples={3}_job={4}{5}".format(scene_file, scene_id,
                                        identifier, samples, job, image_ext)
  else:
    return "scene={0}_{1}_{2}samples={3}{4}".format(scene_file, scene_id,
                                        identifier, samples, image_ext)


def get_scene_filename(scene_file, scene_id, identifier, samples, job,
                       scene_ext = ".pbrt"):
  if job:
    return ".scene={0}_{1}_{2}samples={3}_job={4}{5}".format(scene_file, scene_id,
                                        identifier, samples, job, scene_ext)
  else:
    return ".scene={0}_{1}_{2}samples={3}{4}".format(scene_file, scene_id,
                                        identifier, samples, scene_ext)

def parse_args(args = None):
  parser = argparse.ArgumentParser(description="Run pbrt jobs on antill.")
  
  parser.add_argument("--anthill", metavar="a", type=bool, default=False)
  parser.add_argument("--scenes", metavar="s", type=str, default="[]")
    
  parser.add_argument("--pbrt", metavar="l", type=str, default="../build/pbrt")
  parser.add_argument("--out_dir", metavar="o", type=str, default="./")

  parser.add_argument("--jobs", metavar="j", type=int, default=1)
  parser.add_argument("--samples", metavar="n", type=int, default=1)
  parser.add_argument("--threads", metavar="t", type=int, default=4)

  parser.add_argument("--xres", type=int, default=600)
  parser.add_argument("--yres", type=int, default=800)
  parser.add_argument("--scale", type=float, default=1)
  parser.add_argument("--integrators", type=str, default="[\"nullpath\"]")
  parser.add_argument("--integrator_flags", metavar="f",  
                      type=str, default="{}")
  parser.add_argument("--custom_scene_file", type=bool, default=False)

  args = parser.parse_args(args) if args else parser.parse_args()

  args.scenes = ast.literal_eval(args.scenes)
  args.integrators = ast.literal_eval(args.integrators)
  args.integrator_flags = json.loads(args.integrator_flags)
  return args

def run_job(pbrt, numThreads, pbrt_file, out_dir, save_output = False, seed = 0, anthill = False, echo = False):
  run_pbrt_command = "{0} --seed {1} --nthreads={2} {3}".format(pbrt, seed, numThreads, pbrt_file)
  
  pbrt_filename = os.path.basename(pbrt_file)
  stats_file = pbrt_filename[1:] if pbrt_filename[0] == "." else pbrt_filename
  output_to_txt = " {0}_stats.txt".format(
	os.path.abspath(os.path.join(out_dir, os.path.basename(stats_file)))  )
 
  if echo:
    print(["./run_job.sh", str(run_pbrt_command), str(output_to_txt)])

  if anthill:
    subprocess.call(["./run_anthill_job.sh", 
                    pbrt_file, args.num_threads, 
                    run_pbrt_command, output_to_txt],
                    cwd = ".")
  else:
    if save_output:
      subprocess.call(["./run_job.sh", run_pbrt_command, output_to_txt],
                     cwd = ".")
    else:
      subprocess.call([pbrt, pbrt_file], cwd=".") 

if __name__ == "__main__":
  # Parse arguments, input shown below:
  args = parse_args()
  if len(os.sys.argv) < 2: 
    args.pbrt = "../build/pbrt"
    args.out_dir = "./"
    args.jobs = 1
    args.samples = 32
    args.xres = 600
    args.yres = 800
    args.cropwindow = "0 1.0 0 1.0"
    args.save_output = True
    args.seed = 0
    args.anthill = False

  args_list = []
  args.seed = 0 #1924
  # PAPER SCENES
  args_list = [
    #setup_veach_comparison_one(deepcopy(args)),	# Figure 2: NEE bad (forward scatter, thin)
    #setup_veach_comparison_two(deepcopy(args)),	# Figure 2: Undirectional bad (backward scater, dense)
    #setup_ea_cloud(deepcopy(args)),                    # Figure 3: EA in Cloud
    #setup_colored_smoke_one(deepcopy(args)),           # Figure 4: Colored smoke with only two scattering channels.
    #setup_colored_smoke_one_ms(deepcopy(args)),        # Figure 4: Colored smoke with three scattering channels.
    #setup_colored_smoke_one_new(deepcopy(args)),           # Figure 4: Colored smoke with only two scattering channels.
    #setup_colored_smoke_one_ms_new(deepcopy(args)),        # Figure 4: Colored smoke with three scattering channels.
    #setup_colored_smoke_two(deepcopy(args)),           # Figure 4: Colored smoke with three scattering channels.
    #setup_colored_smoke_two_ms(deepcopy(args)),        # Figure 4: Colored smoke with three scattering channels.
    #setup_colored_smoke_three(deepcopy(args)),         # Figure 4: Colored smoke with three scattering channels.
    #setup_colored_smoke_three_ms(deepcopy(args))      # Figure 4: Colored smoke with three scattering channels.
    #setup_dense_dragon(deepcopy(args)),                # Figure 5: Dense jade dragon.
    #setup_breakfast(deepcopy(args)),                    # Figure 6: Breakfast Scene.
    #setup_colored_smoke_teaser(deepcopy(args)),        # TEASER: Three colored smoke plumes
    setup_cornell(deepcopy(args)),                      # TEASER: Three colored smoke plumes 
  ]
  
  # FURNACE SCENES
  #args = setup_furnace(args)
  #args = setup_ea_furnace(args)
  #args = setup_furnace_env_light(args)
  #args = setup_furnace_path_view(args)
  #args = setup_infinite(args)

  
  # For each scene, create pbrt jobs to run
  for args in args_list:
      for scene in args.scenes:
        args.scene = scene
        if not args.custom_scene_file:
          scene_files = create_scene_files(args) 
        else:
          scene_files = create_custom_scene_files(args)
        for pbrt_file in scene_files:
          run_job(args.pbrt, args.threads, pbrt_file, args.out_dir, 
                  args.save_output, args.seed, args.anthill, echo = True)

