import unittest
from bifrost_salmonella_subspecies_dtartrate import launcher
from bifrost_salmonella_subspecies_dtartrate import salmonella_subspecies_dtartrate

class Test_salmonella_subspecies_dtartrate(unittest.TestCase):
	def test_bifrost_subspecies(self):
		self.assertEqual(salmonella_subspecies_dtartrate.subspecies("19")[0],"enterica")
		self.assertEqual(salmonella_subspecies_dtartrate.subspecies("19*")[0],"NA")

	def test_bifrost_dtartrate(self):
		self.assertEqual(salmonella_subspecies_dtartrate.dtartrate("/bifrost/test_data/samples/SRR2094561_1.fastq.gz",
																	"/bifrost/test_data/samples/SRR2094561_2.fastq.gz")[0],
																	"")
