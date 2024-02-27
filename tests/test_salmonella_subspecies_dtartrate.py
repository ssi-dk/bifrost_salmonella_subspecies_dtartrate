import os
from pathlib import Path
from tempfile import NamedTemporaryFile

import pytest
import yaml
from bifrost_salmonella_subspecies_dtartrate import (
    NearestSubspecies, launcher, rule__dtartrate,
    salmonella_subspecies_dtartrate)

bifrost_install_dir = Path(os.environ['BIFROST_INSTALL_DIR'])
bifrost_config_and_data_path = Path(f"{bifrost_install_dir}/bifrost/test_data")

component_name = "bifrost_salmonella_subspecies_dtartrate"
component_path = bifrost_install_dir/"bifrost/components"/component_name

@pytest.fixture
def config():
	return yaml.safe_load(open(component_path/component_name/"config.yaml",'r'))

class Object:
	def __init__(self):
		pass

@pytest.fixture
def tmpfile():
	yield NamedTemporaryFile("w+")

class Test_salmonella_subspecies_dtartrate:
	test_dir = bifrost_config_and_data_path / "output/test__salmonella_subspecies_dtartrate"
	def test_bifrost_subspecie_clean_st(self, config):
		mlst_db = component_path / config['resources']['mlst_db']
		subspecies_reference = component_path/config['resources']['subspecies_reference']
		subspecies = NearestSubspecies.subspecies_from_st(19, mlst_db, subspecies_reference)[0][1]
		expected_subspecies = "enterica"
		assert subspecies == expected_subspecies

	def silenced_cause_slow_test_bifrost_dtartrate(self, config, tmpfile):
		input=(bifrost_config_and_data_path/"samples/SRR2094561_1.fastq.gz",
			   bifrost_config_and_data_path/"samples/SRR2094561_2.fastq.gz")
		output = Object()
		output._file = tmpfile.name
		params = Object()
		dtartrate_db = component_path / config['resources']['dtartrate_db']
		params.dtartrate_db = dtartrate_db
		log = Object()
		log.err_file = "logfile.txt"
		rule__dtartrate.rule__dtartrate(input, output, params, log)
		tmpfile.seek(0)
		result = tmpfile.read()
		assert result == "{'C': 32}"
