import pytest
from bifrost_salmonella_subspecies_dtartrate import launcher
from bifrost_salmonella_subspecies_dtartrate import salmonella_subspecies_dtartrate
from bifrost_salmonella_subspecies_dtartrate import NearestSubspecies
from pathlib import Path
import yaml
import os

bifrost_install_dir = Path(os.environ['BIFROST_INSTALL_DIR'])
bifrost_config_and_data_path = Path(f"{bifrost_install_dir}/bifrost/test_data")

component_name = "bifrost_salmonella_subspecies_dtartrate"
component_path = bifrost_install_dir/"bifrost/components"/component_name

@pytest.fixture
def config():
	return yaml.safe_load(open(component_path/component_name/"config.yaml",'r'))

class Test_salmonella_subspecies_dtartrate:
	test_dir = bifrost_config_and_data_path / "output/test__salmonella_subspecies_dtartrate"
	def test_bifrost_subspecie_clean_st(self, config):
		mlst_db = component_path / config['resources']['mlst_db']
		print(f"component_path: {component_path}")
		print(f"mlst_db: {mlst_db}")
		subspecies_reference = component_path/config['resources']['subspecies_reference']
		subspecies = NearestSubspecies.subspecies_from_st(19, mlst_db, subspecies_reference)[0][1]
		expected_subspecies = "enterica"
		assert subspecies == expected_subspecies

	# def test_bifrost_dtartrate(self, config):
	# 	dtartrate_db = component_path / config['resources']['dtartrate_db']
	# 	assert salmonella_subspecies_dtartrate.dtartrate(
	# 		bifrost_config_and_data_path/"samples/SRR2094561_1.fastq.gz",
	# 		bifrost_config_and_data_path/"samples/SRR2094561_2.fastq.gz",
	# 		dtartrate_db, component_path)[0] == (
	# 		"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC")
	# 	assert salmonella_subspecies_dtartrate.dtartrate(
	# 		bifrost_config_and_data_path/"samples/2108T25533_S3_L555_R1_001.fastq.gz",
	# 		bifrost_config_and_data_path/"samples/2108T25533_S3_L555_R2_001.fastq.gz",
	# 		dtartrate_db, component_path)[0] == (
	# 		"TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT")
