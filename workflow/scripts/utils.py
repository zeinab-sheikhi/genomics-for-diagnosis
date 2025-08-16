import os

from csv import DictReader
from yaml_validator import load_configfile

config_path = "../config"
config = load_configfile(os.path.join(config_path, "config.yaml"))

