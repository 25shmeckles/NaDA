import pytest
import numpy as np
import douwelib as dl

class TestFastaPrasers:
    
    def setup_method(self):
        file_name = './test_cases/test_fasta_1.fasta'
        self.data = dl.parse_fasta_file(file_name)
        
    def test_data_size(self):
        assert len(self.data) == 5
        
    def test_sequence_ids_are_correct(self):
        for id in '12347':
            assert id in self.data
        
    def test_no_wrong_sequences(self):
        for letter in 'AGTCN':
            assert letter not in self.data.values()
            
class TestScoreStrip:

    def setup_method(self):
        self.path = 'C:/Users/Douwe/Documents/Python/Sequence_files'
        self.data = dl.parse_fasta_file_stripscores(self.path)
    
    def test_data_len(self):
        assert len(self.data) >= 200
        
    def test_no_sequence_and_id(self):
        for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            assert letter not in self.data

class TestQualityScoreConvert:
    
    def setup_method(self):
        sequence_file = 'C:/Users/Douwe/Documents/Python/test_cases/test_fastq2.done_fastq'
        data = dl.parse_fasta_file_error(sequence_file)
        id_ = list(data.keys())[0]
        self.values = dl.convert_qualityscore(data[id_]['score'])
        
    def test_data_is_values(self):
        for character in '!@#$%^&*()?><:;,.[]\=-':
            assert character not in self.values
        
    def test_value_1_or_smaller(self):
        for value in self.values:
            assert value <= 1
        
class TestDoneFastqParser:
    
    def setup_method(self):
        sequence_file = 'C:/Users/Douwe/Documents/Python/test_cases/test_fastq2.done_fastq'
        self.data = dl.parse_fasta_file_error(sequence_file)
        id_ = list(self.data.keys())[0]
        self.score = self.data[id_]['score']

    def check_valid_DNA_sequence(self, s):
        for l in set(s.upper()):
            if not l in 'ACTGN':
                return False
        return True
        
    def test_has_id(self):
        for id in '@':
            assert id in list(self.data.keys())[0]
           
    def test_sequence_correct(self):
        for k, v in self.data.items():
            assert self.check_valid_DNA_sequence(v['sequence']) == True
            
    def test_score_correct(self):
        for letter in 'ABCDEFGHIJKLMNOPQRSTUVWXYZ':
            assert letter not in self.score
            