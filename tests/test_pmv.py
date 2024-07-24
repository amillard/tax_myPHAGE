# Dependencies:
# pip install pytest-mock
import os
import shutil
import pytest
import pandas as pd
from taxmyphage.pmv import ClusteringOnGenomicSimilarity


class TestClusteringOnGenomicSimilarity:

    # Given a file and a reference, when running ClusteringOnGenomicSimilarity, it should create a blastn database, blast against itself, parse the blastn file, calculate distances, cluster all, and return a tuple of a pandas dataframe and a string
    def test_run(self, tmpdir):
        test_file_dir = os.path.dirname(__file__)

        # Create an instance of ClusteringOnGenomicSimilarity
        reference = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "pmv",
                                 "pmv_in.fa"
                                 )
        
        reference_copy = shutil.copyfile(reference, tmpdir / "pmv_in.fa")

        pmv = ClusteringOnGenomicSimilarity(reference_copy, reference, output_dir=tmpdir)

        # Call the run method
        result = pmv.run()

        # Assert the return value
        assert isinstance(result, tuple)
        assert isinstance(result[0], pd.DataFrame)
        assert isinstance(result[1], str)

    # Given a file and a reference, when running ClusteringOnGenomicSimilarity with a genus threshold of 0, it should return a dataframe with all genomes as same cluster
    def test_genus_threshold_zero(self, tmpdir):
        test_file_dir = os.path.dirname(__file__)

        # Create an instance of ClusteringOnGenomicSimilarity
        reference = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "pmv",
                                 "pmv_in.fa"
                                 )
        
        reference_copy = shutil.copyfile(reference, tmpdir / "pmv_in.fa")


        pmv = ClusteringOnGenomicSimilarity(reference_copy, reference, genus_threshold=0, output_dir=tmpdir)

        # Call the run method
        result = pmv.run()

        # Assert the return value
        assert isinstance(result, tuple)
        assert isinstance(result[0], pd.DataFrame)
        assert isinstance(result[1], str)

        # Assert the dataframe contains all genomes as singletons
        df = result[0]
        assert len(df) == len(pmv.size_dict)
        assert df['genus_cluster'].nunique() == 1

    # Given a file and a reference, when running ClusteringOnGenomicSimilarity with a genus threshold of 100, it should return a dataframe with all genomes in a single genus cluster
    def test_genus_threshold_100(self, tmpdir):
        test_file_dir = os.path.dirname(__file__)

        # Create an instance of ClusteringOnGenomicSimilarity
        reference = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "pmv",
                                 "pmv_in.fa"
                                 )
        
        reference_copy = shutil.copyfile(reference, tmpdir / "pmv_in.fa")


        pmv = ClusteringOnGenomicSimilarity(reference_copy, reference, genus_threshold=100, output_dir=tmpdir)

        # Call the run method
        result = pmv.run()

        # Assert the return value
        assert isinstance(result, tuple)
        assert isinstance(result[0], pd.DataFrame)
        assert isinstance(result[1], str)

        # Assert the dataframe contains all genomes as singletons
        df = result[0]
        assert len(df) == len(pmv.size_dict)
        assert df['genus_cluster'].nunique() == 3

    # Given a file and a reference, when running ClusteringOnGenomicSimilarity with a non-existent file, it should raise a FileNotFoundError
    def test_nonexistent_file(self):
        with pytest.raises(FileNotFoundError):
            pmv = ClusteringOnGenomicSimilarity('nonexistent_file', 'reference')
            pmv.run()

    # Given a file and a reference, when running ClusteringOnGenomicSimilarity with a non-existent reference, it should raise a FileNotFoundError
    def test_nonexistent_reference(self):
        test_file_dir = os.path.dirname(__file__)

        # Create an instance of ClusteringOnGenomicSimilarity
        first_file = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "pmv",
                                 "pmv_in.fa"
                                 )
        
        # Assert that running the ClusteringOnGenomicSimilarity raises a FileNotFoundError
        with pytest.raises(FileNotFoundError):
            pmv = ClusteringOnGenomicSimilarity(first_file, 'nonexistent_reference')
            pmv.run()

    def test_good_df_returned(self, tmpdir):
        test_file_dir = os.path.dirname(__file__)

        # Create an instance of ClusteringOnGenomicSimilarity
        reference = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "pmv",
                                 "pmv_in.fa"
                                 )
        
        reference_copy = shutil.copyfile(reference, tmpdir / "pmv_in.fa")

        pmv = ClusteringOnGenomicSimilarity(reference_copy, reference, output_dir=tmpdir)

        # Call the run method
        result = pmv.run()

        # Assert the return value
        assert isinstance(result, tuple)
        assert isinstance(result[0], pd.DataFrame)
        assert isinstance(result[1], str)

        # Assert the dataframe contains all genomes as singletons
        df = result[0]
        assert len(df) == len(pmv.size_dict)
        assert df['genus_cluster'].nunique() == 1
        assert df['species_cluster'].nunique() == 2

        table_final = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "pmv",
                                 "pmv_in.fa.genus_species_clusters.tsv"
                                 )
    
        table_final = pd.read_csv(table_final, sep='\t')
        assert table_final.equals(df)

####################################################################################################