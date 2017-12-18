import pytest
import douwelib as dl

class TestStripVarianceData:
    
    def setup_method(self):
        file_name = './test_cases/test_vcf.vcf'
        self.data = dl.data_vcf_file(file_name)
        self.id_ = list(self.data.keys())[0]
        self.variance = self.data[self.id_]['variance']
        self.backbone = self.data[self.id_]['backbone']
        
    def test_id_correct_line(self):
        test = True
        if '##' in self.id_:
            test = False
        elif './.' in self.id_:
            test = False
        elif 'GT:AD:DP:GQ' in self.id_:
            test = False
        assert test == True
    
    def test_data_no_empty_sequences(self):
        assert './.' not in self.data
        
    def test_variance_no_empty_sequence(self):
        assert './.' not in self.variance 
        
    def test_backbone_no_empty_sequence(self):
        assert './.' not in self.backbone
        
    def test_variance_has_no_backbone(self):
        assert 'BB' not in self.variance
        
class TestStripWholeSequence:
    
    def setup_method(self):
        file_name = './test_cases/test_vcf.vcf'
        self.data = dl.vcf_whole_sequence_strip(file_name)      
    
    def test_only_sequence(self):
        assert '#' not in self.data
        
    def test_does_contain_seq(self):
        assert './.' in self.data[0]

class TestMutatedVarianceOnly:
    
    def setup_method(self):
        file_name = './test_cases/test_vcf.vcf'
        self.data = dl.data_vcf_file(file_name)
        self.id_ = list(self.data.keys())[0]
        self.variance = self.data[self.id_]['variance']
        self.backbone = self.data[self.id_]['backbone']
        self.data_all = dl.vcf_whole_sequence_strip(file_name)
        self.extended = dl.mutated_reads_vcf_only(self.variance, self.data_all)[0]
        self.highmutated = dl.mutated_reads_vcf_only(self.variance, self.data_all)[1]
   
    def test_len_extended(self):
        for item in self.extended:
            assert len(item) == 5
        
    def test_highmutated(self):
        test = True
        for item in self.highmutated:
            i = item.split('\t', 10)[9]
            items = i.split(':', 4)[1]
            if len(item) > 1:
                test = False
            elif len(items) == 3:
                score = (int(items.split(',', 2)[1]))/(int(i.split(':', 4)[2]))
                if score <= 0.25:
                    test = False                     
        assert test == True