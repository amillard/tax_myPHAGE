import os
import pytest
import pandas as pd

from taxmyphage.classify import new_genus
from taxmyphage.classify import current_genus_current_species
from taxmyphage.classify import current_genus_new_species
from taxmyphage.classify import new_genus_new_species
from taxmyphage.classify import assess_taxonomic_info
from taxmyphage.classify import classification
from taxmyphage.classify import classification_viridic

from taxmyphage.handle_files import read_VMR


####################################################################################################

class TestNewGenus:
    """
    Test the new_genus function
    """

    # Returns a dictionary with taxonomic information when given a valid query genus cluster number and summary output path
    def test_valid_query_genus_cluster_number(self):
        """
        Returns a dictionary with taxonomic information when given a valid query genus cluster number and summary output path
        """
        # Arrange
        query_genus_cluster_number = 0
        dict_genus_cluster_2_genus_name = {1: "Genus1"}
        summary_output_path = ""
        message = "Query does not fall within a current genus or species as defined by ICTV"
    
        # Act
        result = new_genus(query_genus_cluster_number, dict_genus_cluster_2_genus_name, summary_output_path, message)
    
        # Assert
        assert isinstance(result, dict)
        assert result["Genus"] == "New_genus"
        assert result["Species"] == "New_species"

    # Prints a result message when given an invalid query genus cluster number and summary output path
    def test_empty_dict_genus_cluster_genus_name(self):
        """
        Prints a result message when given an invalid query genus cluster number and summary output path
        """
        # Arrange
        query_genus_cluster_number = 0
        dict_genus_cluster_2_genus_name = {}
        summary_output_path = ""
        message = "Query does not fall within a current genus or species as defined by ICTV"
    
        # Act
        result = new_genus(query_genus_cluster_number, dict_genus_cluster_2_genus_name, summary_output_path, message)
    
        # Assert
        assert isinstance(result, dict)
        assert result["Genus"] == "New_genus"
        assert result["Species"] == "New_species"

    def test_writes_taxonomic_information_to_summary_output_file_if_provided(self, tmp_path):
        """
        writes the taxonomic information to a summary output file if provided
        """

        # Arrange
        summary_output_path = tmp_path / "summary.txt"
        query_genus_cluster_number = 0
        dict_genus_cluster_2_genus_name = {}
        message = "Query does not fall within a current genus or species as defined by ICTV"

        # Act
        new_genus(query_genus_cluster_number, dict_genus_cluster_2_genus_name, summary_output_path, message)
    

        # Assert
        assert summary_output_path.exists()
        with open(summary_output_path, "r", encoding="utf-8") as file:
            content = file.read()
            assert "Genus: New genus" in content
            assert "Species: New species" in content
            assert f"INFO: {message}" in content

####################################################################################################

class TestCurrentGenusCurrentSpecies:
    """
    Test the current_genus_current_species function
    """

    # correctly identifies the predicted genus and species
    def test_correctly_identifies_predicted_genus_and_species(self):
        """
        correctly identifies the predicted genus and species
        """
        
        # Arrange
        query_species_cluster_number = 1
        dict_species_cluster_2_species_name = {1: "species1"}
        summary_output_path = ""
        dict_genus_cluster_2_genus_name = {1: "genus1"}
        query_genus_cluster_number = 1
        merged_df = pd.DataFrame({"Species": ["species1"], "Genus": ["genus1"], "Class": ["class1"], "Family": ["family1"], "Subfamily": ["subfamily1"]})
        mash_df = pd.DataFrame()

        message = "Current ICTV taxonomy and VIRIDIC-algorithm output appear to be consistent at the genus level"

        # Act
        result = current_genus_current_species(
            query_species_cluster_number,
            dict_species_cluster_2_species_name,
            summary_output_path,
            dict_genus_cluster_2_genus_name,
            query_genus_cluster_number,
            merged_df,
            mash_df,
            message,
        )

        # Assert
        assert result["Genus"] == "genus1"
        assert result["Species"] == "species1"
        assert result["Message"] == message

    # prints the taxonomic information of the matching species row
    def test_prints_taxonomic_information_of_matching_species_row(self, capsys):
        """
        prints the taxonomic information of the matching species row
        """

        # Arrange
        query_species_cluster_number = 1
        dict_species_cluster_2_species_name = {1: "species1"}
        summary_output_path = ""
        dict_genus_cluster_2_genus_name = {1: "genus1"}
        query_genus_cluster_number = 1
        merged_df = pd.DataFrame({"Species": ["species1"], "Genus": ["genus1"], "Class": ["class1"], "Family": ["family1"], "Subfamily": ["subfamily1"]})
        mash_df = pd.DataFrame()

        message = "Current ICTV taxonomy and VIRIDIC-algorithm output appear to be consistent at the genus level"

        # Act
        current_genus_current_species(
            query_species_cluster_number,
            dict_species_cluster_2_species_name,
            summary_output_path,
            dict_genus_cluster_2_genus_name,
            query_genus_cluster_number,
            merged_df,
            mash_df,
            message,
        )

        # Assert
        captured = capsys.readouterr()
        assert "QUERY is in the genus: genus1 and is species: species1" in captured.out
        assert "Class: class1" in captured.out
        assert "Family: family1" in captured.out
        assert "Subfamily: subfamily1" in captured.out
        assert "Genus: genus1" in captured.out
        assert "Species: species1" in captured.out

    # writes the taxonomic information to a summary output file if provided
    def test_writes_taxonomic_information_to_summary_output_file_if_provided(self, tmp_path):
        """
        writes the taxonomic information to a summary output file if provided
        """

        # Arrange
        query_species_cluster_number = 1
        dict_species_cluster_2_species_name = {1: "species1"}
        summary_output_path = tmp_path / "summary.txt"
        dict_genus_cluster_2_genus_name = {1: "genus1"}
        query_genus_cluster_number = 1
        merged_df = pd.DataFrame({"Species": ["species1"], "Genus": ["genus1"], "Class": ["class1"], "Family": ["family1"], "Subfamily": ["subfamily1"]})
        mash_df = pd.DataFrame()

        message = "Current ICTV taxonomy and VIRIDIC-algorithm output appear to be consistent at the genus level"

        # Act
        current_genus_current_species(
            query_species_cluster_number,
            dict_species_cluster_2_species_name,
            summary_output_path,
            dict_genus_cluster_2_genus_name,
            query_genus_cluster_number,
            merged_df,
            mash_df,
            message,
        )

        # Assert
        assert summary_output_path.exists()
        with open(summary_output_path, "r", encoding="utf-8") as file:
            content = file.read()
            assert "Class: class1" in content
            assert "Family: family1" in content
            assert "Subfamily: subfamily1" in content
            assert "Genus: genus1" in content
            assert "Species: species1" in content
            assert f"INFO: {message}" in content

    # handles empty input values for query_species_cluster_number and query_genus_cluster_number
    def test_return_correct_directory(self):
        """
        handles empty input values for query_species_cluster_number and query_genus_cluster_number
        """

        # Arrange
        query_species_cluster_number = 1
        dict_species_cluster_2_species_name = {1: "species1"}
        summary_output_path = ""
        dict_genus_cluster_2_genus_name = {1: "genus1"}
        query_genus_cluster_number = 1
        merged_df = pd.DataFrame({"Species": ["species1"], "Genus": ["genus1"], "Class": ["class1"], "Family": ["family1"], "Subfamily": ["subfamily1"]})
        mash_df = pd.DataFrame()

        message = "Current ICTV taxonomy and VIRIDIC-algorithm output appear to be consistent at the genus level"

        # Act
        result = current_genus_current_species(
            query_species_cluster_number,
            dict_species_cluster_2_species_name,
            summary_output_path,
            dict_genus_cluster_2_genus_name,
            query_genus_cluster_number,
            merged_df,
            mash_df,
            message,
        )

        # Assert
        assert result == {
            "Class": "class1",
            "Family": "family1",
            "Subfamily": "subfamily1",
            "Genus": "genus1",
            "Species": "species1",
            "Message": message,
        }

####################################################################################################

class TestCurrentGenusNewSpecies:
    """
    Test the current_genus_new_species function
    """

    # Returns a dictionary containing taxonomic information for a query genome that is within a current genus but represents a new species
    def test_returns_dictionary_with_taxonomic_information(self):
        """
        Returns a dictionary containing taxonomic information for a query genome that is within a current genus but represents a new species
        """

        # Arrange
        summary_output_path = ""
        dict_genus_cluster_2_genus_name = {1: "Genus1", 2: "Genus2"}
        query_genus_cluster_number = 1
        merged_df = pd.DataFrame({"Genus": ["Genus1"], "Class": ["Class1"], "Family": ["Family1"], "Subfamily": ["Subfamily1"]})
        mash_df = pd.DataFrame()

        message = "The number of expected genera is different from the predicted number of genus clusters. It will require more manual curation"

        # Act
        result = current_genus_new_species(summary_output_path,
                                           dict_genus_cluster_2_genus_name,
                                           query_genus_cluster_number,
                                           merged_df,
                                           mash_df,
                                           message,
                                           )

        # Assert
        expected_result = {
            "Class": "Class1",
            "Family": "Family1",
            "Subfamily": "Subfamily1",
            "Genus": "Genus1",
            "Species": "Genus1 new_name",
            "Message": message,
        }
        assert result == expected_result

    # Prints the predicted genus and the taxonomic information for the query genome
    def test_prints_predicted_genus_and_taxonomic_information(self, capsys):
        """
        Prints the predicted genus and the taxonomic information for the query genome
        """

        # Arrange
        summary_output_path = ""
        dict_genus_cluster_2_genus_name = {1: "Genus1", 2: "Genus2"}
        query_genus_cluster_number = 1
        merged_df = pd.DataFrame({"Genus": ["Genus1"], "Class": ["Class1"], "Family": ["Family1"], "Subfamily": ["Subfamily1"]})
        mash_df = pd.DataFrame()

        message = "The number of expected genera is different from the predicted number of genus clusters. It will require more manual curation"

        # Act
        current_genus_new_species(summary_output_path, 
                                  dict_genus_cluster_2_genus_name, 
                                  query_genus_cluster_number, 
                                  merged_df,
                                  mash_df,
                                  message,
                                )

        # Assert
        captured = capsys.readouterr()
        assert "Query sequence is:" in captured.out
        assert "Class: Class1" in captured.out
        assert "Family: Family1" in captured.out
        assert "Subfamily: Subfamily1" in captured.out
        assert "Genus: Genus1" in captured.out
        assert "Species: Genus1 new_name" in captured.out

    # Writes the taxonomic information to a summary output file if the path is provided
    def test_writes_taxonomic_information_to_summary_output_file(self, tmp_path):
        """
        Writes the taxonomic information to a summary output file if the path is provided
        """
        
        # Arrange
        summary_output_path = tmp_path / "summary_output.txt"
        dict_genus_cluster_2_genus_name = {1: "Genus1", 2: "Genus2"}
        query_genus_cluster_number = 1
        merged_df = pd.DataFrame({"Genus": ["Genus1"], "Class": ["Class1"], "Family": ["Family1"], "Subfamily": ["Subfamily1"]})
        mash_df = pd.DataFrame()

        message = "The number of expected genera is different from the predicted number of genus clusters. It will require more manual curation"

        # Act
        current_genus_new_species(summary_output_path,
                                  dict_genus_cluster_2_genus_name,
                                  query_genus_cluster_number,
                                  merged_df,
                                  mash_df,
                                  message,
                                  )

        # Assert
        assert summary_output_path.exists()
        with open(summary_output_path, "r", encoding="utf-8") as file:
            content = file.read()
            assert "Class: Class1" in content
            assert "Family: Family1" in content
            assert "Subfamily: Subfamily1" in content
            assert "Genus: Genus1" in content
            assert "Species: Genus1 new_name" in content
            assert f"INFO: {message}" in content


    # Raises a KeyError if the query genus cluster number is not found in the dictionary of genus cluster to genus name
    def test_raises_key_error_if_query_genus_cluster_number_not_found_in_dictionary(self):
        """
        Raises a KeyError if the query genus cluster number is not found in the dictionary of genus cluster to genus name
        """
        # Arrange
        summary_output_path = ""
        dict_genus_cluster_2_genus_name = {1: "Genus1", 2: "Genus2"}
        query_genus_cluster_number = 3
        merged_df = pd.DataFrame({"Genus": ["Genus1"], "Class": ["Class1"], "Family": ["Family1"], "Subfamily": ["Subfamily1"]})
        mash_df = pd.DataFrame()

        # Act and Assert
        with pytest.raises(KeyError):
            current_genus_new_species(summary_output_path,
                                      dict_genus_cluster_2_genus_name,
                                      query_genus_cluster_number,
                                      merged_df,
                                      mash_df,
                                      message='',
                                      )

####################################################################################################
class TestNewGenusNewSpecies:
    """
    Test the new_genus_new_species function
    """

    # Returns a dictionary of taxonomic information with 'New_genus' and 'New_species' as the genus and species values respectively.
    def test_returns_dictionary_with_new_genus_and_species(self):
        """
        Returns a dictionary of taxonomic information with 'New_genus' and 'New_species' as the genus and species values respectively.
        """
        # Arrange
        summary_output_path = ""
        mash_df = pd.DataFrame()
    
        message = "Query is a new genus and species. You could try running again with if you larger distance"

        # Act
        result = new_genus_new_species(summary_output_path,
                                       mash_df,
                                       message
                                       )
    
        # Assert
        assert isinstance(result, dict)
        assert result["Genus"] == "New_genus"
        assert result["Species"] == "New_species"
        assert result["Message"] == message

    # Prints a message indicating that the query sequence is likely the first representative of both a new species and new genus.
    def test_prints_message_indicating_new_genus_and_species(self, capsys):
        """
        Prints a message indicating that the query sequence is likely the first representative of both a new species and new genus.
        """

        # Arrange
        summary_output_path = ""
        mash_df = pd.DataFrame()
    
        message = "Query is a new genus and species. You could try running again with if you larger distance"

        # Act
        new_genus_new_species(summary_output_path, mash_df, message)
        captured = capsys.readouterr()
    
        # Assert
        assert "Query does not fall within a current genus or species as defined by ICTV" in captured.out
        assert "Therefore the query sequence is likely the first representative of both a new species and new genus." in captured.out

    # Writes a summary output file with a message indicating that the query sequence can not be classified within a current genus or species.
    def test_writes_summary_output_file_with_message(self, tmp_path):
        """
        Writes a summary output file with a message indicating that the query sequence can not be classified within a current genus or species.
        """

        # Arrange
        summary_output_path = tmp_path / "summary_output.txt"
        mash_df = pd.DataFrame()

        message = "Query is a new genus and species. You could try running again with if you larger distance"

        # Act
        new_genus_new_species(summary_output_path, mash_df, message)

        # Assert
        assert os.path.exists(summary_output_path)

        with open(summary_output_path, "r", encoding="utf-8") as file:
            content = file.read()
            assert "Query sequence can not be classified within a current genus or species" in content

    # Returns a dictionary of taxonomic information with 'Unknown' as the values for all taxonomic ranks except for 'Genus' and 'Species'.
    def test_returns_dictionary_with_unknown_taxonomic_ranks(self):
        """
        Returns a dictionary of taxonomic information with 'Unknown' as the values for all taxonomic ranks except for 'Genus' and 'Species'.
        """

        # Arrange
        summary_output_path = ""
        mash_df = pd.DataFrame()
    
        message = "Query is a new genus and species. You could try running again with if you larger distance"

        # Act
        result = new_genus_new_species(summary_output_path, mash_df, message)
    
        # Assert
        assert isinstance(result, dict)
        assert result["Realm"] == "Unknown"
        assert result["Kingdom"] == "Unknown"
        assert result["Phylum"] == "Unknown"
        assert result["Class"] == "Unknown"
        assert result["Order"] == "Unknown"
        assert result["Family"] == "Unknown"
        assert result["Subfamily"] == "Unknown"
        assert result["Message"] == message

    # Does not compare against all other known phages, only those that have been classified.
    def test_does_not_compare_against_all_known_phages(self, capsys):
        """
        Does not compare against all other known phages, only those that have been classified.
        """

        # Arrange
        summary_output_path = ""
        mash_df = pd.DataFrame()
    
        # Act
        new_genus_new_species(summary_output_path, mash_df, message="")
        captured = capsys.readouterr()
    
        # Assert
        assert "WARNING:: taxmyphage does not compare against all other known phages" in captured.out

####################################################################################################

class TestAssessTaxonomicInfo:
    """
    Test the assess_taxonomic_info function
    """

    # Classifies query genome as current genus and current species
    def test_current_genus_current_species_with_class_key(self):
        """
        Classifies query genome as current genus and current species
        """

        query_genus_cluster_number = 1
        query_species_cluster_number = 2
        list_ICTV_genus_clusters = [1, 3, 5]
        list_ICTV_species_clusters = [2, 4, 6]
        dict_genus_cluster_2_genus_name = {1: "Genus1", 2: "Genus2", 3: "Genus3"}
        dict_species_cluster_2_species_name = {1: "Species1", 2: "Species2", 3: "Species3"}
        summary_output_path = ""
        merged_df = pd.DataFrame({"Species": ["Species1", "Species2", "Species3"], "Class": ["", "", ""], "Family": ["", "", ""], "Subfamily": ["", "", ""], "Genus": ["Genus1", "Genus2", "Genus3"]})
        mash_df = pd.DataFrame()

        message = "Current ICTV taxonomy and VIRIDIC-algorithm output appear to be consistent at the genus level"

        taxonomic_info = assess_taxonomic_info(
            query_genus_cluster_number=query_genus_cluster_number,
            query_species_cluster_number=query_species_cluster_number,
            list_ICTV_genus_clusters=list_ICTV_genus_clusters,
            list_ICTV_species_clusters=list_ICTV_species_clusters,
            dict_genus_cluster_2_genus_name=dict_genus_cluster_2_genus_name,
            dict_species_cluster_2_species_name=dict_species_cluster_2_species_name,
            summary_output_path=summary_output_path,
            merged_df=merged_df,
            mash_df=mash_df,
            message=message,
        )

        assert taxonomic_info == {
            "Class": "",
            "Family": "",
            "Subfamily": "",
            "Genus": "Genus2",
            "Species": "Species2",
            "Message": message,
        }

    # Classifies query genome as current genus and new species with fixed merged_df
    def test_current_genus_new_species_with_merged_df_row_fixed(self):
        """
        Classifies query genome as current genus and new species with fixed merged_df
        """

        query_genus_cluster_number = 1
        query_species_cluster_number = 4
        list_ICTV_genus_clusters = [1, 3, 5]
        list_ICTV_species_clusters = [2, 4, 6]
        dict_genus_cluster_2_genus_name = {1: "Genus1", 2: "Genus2", 3: "Genus3"}
        dict_species_cluster_2_species_name = {1: "Species1", 2: "Species2", 3: "Species3", 4: "Species4"}
        summary_output_path = ""
        merged_df = pd.DataFrame({"Species": ["Species4"], "Class": ["Class1"], "Family": ["Family1"], "Subfamily": ["Subfamily1"], "Genus": ["Genus1"]})
        mash_df = pd.DataFrame()

        message = "The number of expected genera is different from the predicted number of genus clusters. It will require more manual curation"

        taxonomic_info = assess_taxonomic_info(
            query_genus_cluster_number=query_genus_cluster_number,
            query_species_cluster_number=query_species_cluster_number,
            list_ICTV_genus_clusters=list_ICTV_genus_clusters,
            list_ICTV_species_clusters=list_ICTV_species_clusters,
            dict_genus_cluster_2_genus_name=dict_genus_cluster_2_genus_name,
            dict_species_cluster_2_species_name=dict_species_cluster_2_species_name,
            summary_output_path=summary_output_path,
            merged_df=merged_df,
            mash_df=mash_df,
            message=message,
        )

        assert taxonomic_info == {
            "Class": "Class1",
            "Family": "Family1",
            "Subfamily": "Subfamily1",
            "Genus": "Genus1",
            "Species": "Species4",
            "Message": message,
        }

    # Query genus cluster number not in list_ICTV_genus_clusters but query species cluster number in list_ICTV_species_clusters
    def test_query_genus_cluster_not_in_list_ICTV_genus_clusters(self):
        """
        Query genus cluster number not in list_ICTV_genus_clusters but query species cluster number in list_ICTV_species_clusters
        """

        query_genus_cluster_number = 1
        query_species_cluster_number = 2
        list_ICTV_genus_clusters = [3, 5]
        list_ICTV_species_clusters = [2, 4, 6]
        dict_genus_cluster_2_genus_name = {1: "Genus1", 2: "Genus2", 3: "Genus3"}
        dict_species_cluster_2_species_name = {1: "Species1", 2: "Species2", 3: "Species3"}
        summary_output_path = ""
        merged_df = pd.DataFrame({"Species": ["Species1", "Species2", "Species3"], "Class": ["", "", ""], "Family": ["", "", ""], "Subfamily": ["", "", ""], "Genus": ["Genus1", "Genus2", "Genus3"]})
        mash_df = pd.DataFrame()

        message = "The number of expected genera is different from the predicted number of genus clusters. It will require more manual curation"

        with pytest.raises(UnboundLocalError):
            assess_taxonomic_info(
                query_genus_cluster_number=query_genus_cluster_number,
                query_species_cluster_number=query_species_cluster_number,
                list_ICTV_genus_clusters=list_ICTV_genus_clusters,
                list_ICTV_species_clusters=list_ICTV_species_clusters,
                dict_genus_cluster_2_genus_name=dict_genus_cluster_2_genus_name,
                dict_species_cluster_2_species_name=dict_species_cluster_2_species_name,
                summary_output_path=summary_output_path,
                merged_df=merged_df,
                mash_df=mash_df,
                message=message,
            )

    # Query genus cluster number not in list_ICTV_genus_clusters and query species cluster number not in list_ICTV_species_clusters
    def test_query_genus_cluster_number_not_in_list_ICTV_genus_clusters_and_query_species_cluster_number_not_in_list_ICTV_species_clusters(self):
        """
        Query genus cluster number not in list_ICTV_genus_clusters and query species cluster number not in list_ICTV_species_clusters
        """
        
        query_genus_cluster_number = 1
        query_species_cluster_number = 2
        list_ICTV_genus_clusters = [3, 4, 5]
        list_ICTV_species_clusters = [6, 7, 8]
        dict_genus_cluster_2_genus_name = {1: "Genus1", 2: "Genus2", 3: "Genus3"}
        dict_species_cluster_2_species_name = {1: "Species1", 2: "Species2", 3: "Species3"}
        summary_output_path = ""
        merged_df = pd.DataFrame({"Species": ["Species1", "Species2", "Species3"], "Class": ["", "", ""], "Family": ["", "", ""], "Subfamily": ["", "", ""], "Genus": ["Genus1", "Genus2", "Genus3"]})
        mash_df = pd.DataFrame()

        message = "Query is a new genus and species. You could try running again with if you larger distance"

        taxonomic_info = assess_taxonomic_info(
            query_genus_cluster_number=query_genus_cluster_number,
            query_species_cluster_number=query_species_cluster_number,
            list_ICTV_genus_clusters=list_ICTV_genus_clusters,
            list_ICTV_species_clusters=list_ICTV_species_clusters,
            dict_genus_cluster_2_genus_name=dict_genus_cluster_2_genus_name,
            dict_species_cluster_2_species_name=dict_species_cluster_2_species_name,
            summary_output_path=summary_output_path,
            merged_df=merged_df,
            mash_df=mash_df,
            message=message,
        )

        assert taxonomic_info == {
            "Realm": "Unknown",
            "Kingdom": "Unknown",
            "Phylum": "Unknown",
            "Class": "Unknown",
            "Order": "Unknown",
            "Family": "Unknown",
            "Subfamily": "Unknown",
            "Genus": "New_genus",
            "Species": "New_species",
            "Message": message,
        }

    # Query genus cluster number in list_ICTV_genus_clusters but query species cluster number not in list_ICTV_species_clusters
    def test_query_genus_cluster_number_in_list_ICTV_genus_clusters_but_query_species_cluster_number_not_in_list_ICTV_species_clusters(self):
        """
        Query genus cluster number in list_ICTV_genus_clusters but query species cluster number not in list_ICTV_species_clusters
        """

        query_genus_cluster_number = 1
        query_species_cluster_number = 2
        list_ICTV_genus_clusters = [1, 3, 5]
        list_ICTV_species_clusters = [2, 4, 6]
        dict_genus_cluster_2_genus_name = {1: "Genus1", 2: "Genus2", 3: "Genus3"}
        dict_species_cluster_2_species_name = {1: "Species1", 2: "Species2", 3: "Species3"}
        summary_output_path = ""
        merged_df = pd.DataFrame({"Species": ["Species1", "Species2", "Species3"], "Class": ["", "", ""], "Family": ["", "", ""], "Subfamily": ["", "", ""], "Genus": ["Genus1", "Genus2", "Genus3"]})
        mash_df = pd.DataFrame()

        message = "The number of expected genera is different from the predicted number of genus clusters. It will require more manual curation"

        taxonomic_info = assess_taxonomic_info(
            query_genus_cluster_number=query_genus_cluster_number,
            query_species_cluster_number=query_species_cluster_number,
            list_ICTV_genus_clusters=list_ICTV_genus_clusters,
            list_ICTV_species_clusters=list_ICTV_species_clusters,
            dict_genus_cluster_2_genus_name=dict_genus_cluster_2_genus_name,
            dict_species_cluster_2_species_name=dict_species_cluster_2_species_name,
            summary_output_path=summary_output_path,
            merged_df=merged_df,
            mash_df=mash_df,
            message=message,
        )

        assert taxonomic_info == {
            "Class": "",
            "Family": "",
            "Subfamily": "",
            "Genus": "Genus2",
            "Species": "Species2",
            "Message": message,
        }
 

####################################################################################################
class TestClassification:
    """
    Test the classification function
    """
    # Classifies the query genome with consistent ICTV and VIRIDIC-algorithm predictions at the genus level
    def test_consistent_predictions_with_genome_column(self, tmp_path):
        """
        Classifies the query genome with consistent ICTV and VIRIDIC-algorithm predictions at the genus level
        """

        # Initialize the required input dataframes and variables
        merged_df = pd.DataFrame({'genus_cluster': [1, 2, 3],
                                  'species_cluster': [1, 2, 3], 
                                  'Genus': ['Genus1', 'Genus2', 'Genus3'], 
                                  'Species': ['Species1', 'Species2', 'Species3'], 
                                  "Class": ["Class1", "Class2", "Class3"], 
                                  "Order": ["Order1", "Order2", "Order3"], 
                                  "Phylum": ["Phylum1", "Phylum2", "Phylum3"], 
                                  "Family": ["Family1", "Family2", "Family3"], 
                                  "Subfamily": ["Subfamily1", "Subfamily2", "Subfamily3"],
                                  'genome': ['genome1.fasta', 'genome2.fasta', 'genome3.fasta'], 
                                  'closest_genome': ['genome1.fasta', 'genome2.fasta', 'genome3.fasta']
                                  })
        query_merged_df = pd.DataFrame({'genus_cluster': [1], 'Genus': ['Genus1'], 'species_cluster': [1], 'genome': ['query_genome.fasta']})
        results_path = tmp_path
        mash_df = pd.DataFrame({'Reference': ['genome1.fasta'],
                                'Query': ['query_genome.fasta'],
                                'distance': [0.197292],
                                'p-value': [4.879559999999999e-267],
                                'shared-hashes':['40/5000'],
                                'ANI':['']
                                })
        prefix = "output_"
        closest_genome = "genome1.fasta"

        # Invoke the classification function
        taxonomic_info = classification(
            merged_df,
            query_merged_df,
            results_path,
            mash_df,
            prefix,
            closest_genome
        )

        print(taxonomic_info)

        # Assert that the taxonomic information is as expected
        assert taxonomic_info["Message"] == "Current ICTV taxonomy and VIRIDIC-algorithm output appear to be consistent at the genus level"
        assert taxonomic_info["Genus"] == "Genus1"
        assert taxonomic_info["Species"] == "Species1"
        assert taxonomic_info["Family"] == "Family1"
        assert taxonomic_info["Order"] == "Order1"
        assert taxonomic_info["Class"] == "Class1"
        assert taxonomic_info["Phylum"] == "Phylum1"

    # Classifies the query genome with inconsistent ICTV and VIRIDIC-algorithm predictions at the genus level
    def test_inconsistent_predictions_with_genus_cluster_column(self, tmp_path):
        """
        Classifies the query genome with inconsistent ICTV and VIRIDIC-algorithm predictions at the genus level
        """
        
        # Initialize the required input dataframes and variables
        merged_df = pd.DataFrame({'Genus': ['Genus1', 'Genus1', 'Genus1'], 
                                  'Species': ['Species1', 'Species1', 'Species1'], 
                                  "Class": ["Class1", "Class2", "Class3"], 
                                  "Order": ["Order1", "Order2", "Order3"], 
                                  "Phylum": ["Phylum1", "Phylum2", "Phylum3"], 
                                  "Family": ["Family1", "Family2", "Family3"], 
                                  "Subfamily": ["Subfamily1", "Subfamily2", "Subfamily3"],
                                  'closest_genome': ['genome1', 'genome1', 'genome1'],
                                  'genus_cluster': [1, 2, 3], 
                                  'species_cluster': [1, 2, 3], 
                                  'genome': ['genome1', 'genome2', 'genome3']
                                  })
        query_merged_df = pd.DataFrame({'genus_cluster': [1], 'Genus': ['Genus1'], 'species_cluster': [1], 'genome': ['query_seq.fasta']})

        results_path = tmp_path
        mash_df = pd.DataFrame({'Reference': ['genome1.fasta'],
                                'Query': ['query_seq.fasta'],
                                'distance': [0.197292],
                                'p-value': [4.879559999999999e-267],
                                'shared-hashes':['40/5000'],
                                'ANI':['']
                                })
        prefix = "output_"
        closest_genome = "genome1"

        # Invoke the classification function
        taxonomic_info = classification(
            merged_df,
            query_merged_df,
            results_path,
            mash_df,
            prefix,
            closest_genome
        )

        message = "The number of expected genera is different from the predicted number of genus clusters. It will require more manual curation"

        # Assert that the taxonomic information is as expected
        assert taxonomic_info["Message"] == message
        assert taxonomic_info["Genus"] == "Genus1"
        assert taxonomic_info["Species"] == "Species1"
        assert taxonomic_info["Family"] == "Family1"
        assert taxonomic_info["Order"] == "Order1"
        assert taxonomic_info["Class"] == "Class1"
        assert taxonomic_info["Phylum"] == "Phylum1"

    # Classifies the query genome as a new genus when the 'genus_cluster' column is missing in the input dataframe
    def test_new_genus_with_missing_genus_cluster_column(self, tmp_path):
        """
        Classifies the query genome as a new genus when the 'genus_cluster' column is missing in the input dataframe
        """

        # Initialize the required input dataframes and variables
        merged_df = pd.DataFrame({'Genus': ['Genus1', 'Genus1', 'Genus1'], 
                                  'Species': ['Species1', 'Species1', 'Species1'], 
                                  "Class": ["Class1", "Class2", "Class3"], 
                                  "Order": ["Order1", "Order2", "Order3"], 
                                  "Phylum": ["Phylum1", "Phylum2", "Phylum3"], 
                                  "Family": ["Family1", "Family2", "Family3"], 
                                  "Subfamily": ["Subfamily1", "Subfamily2", "Subfamily3"],
                                  'closest_genome': ['genome1', 'genome1', 'genome1'],
                                  'species_cluster': [1, 2, 3], 
                                  'genome': ['genome1', 'genome2', 'genome3']
                                  })
        query_merged_df = pd.DataFrame({'genus_cluster': [1], 'Genus': ['Genus1'], 'species_cluster': [1], 'genome': ['query_seq.fasta']})

        results_path = tmp_path
        mash_df = pd.DataFrame({'Reference': ['genome1.fasta'],
                                'Query': ['query_seq.fasta'],
                                'distance': [0.197292],
                                'p-value': [4.879559999999999e-267],
                                'shared-hashes':['40/5000'],
                                'ANI':['']
                                })
        
        prefix = "output_"
        closest_genome = "closest_genome.fasta"

        # Invoke the classification function
        with pytest.raises(KeyError):
            classification(
                merged_df,
                query_merged_df,
                results_path,
                mash_df,
                prefix,
                closest_genome
            )


####################################################################################################

class TestClassificationViridic:

    # Runs VIRIDIC on input fasta files and generates a merged dataframe
    def test_runs_viridic_and_generates_merged_dataframe(self, tmpdir):
        test_file_dir = os.path.dirname(__file__)
        
        # Initialize input variables
        known_taxa_path = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "known_taxa.fa"
                                 )
        query = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "query.fasta"
                                 )
        
        taxa_df = read_VMR(
            os.path.join(test_file_dir,
                        "dummy_genome",
                        "VMR.dummy.xlsx"
                    )
            )
        
        taxa_df = taxa_df.rename(columns={"Accession": "Genbank"})

        taxa_csv_output_path = os.path.join(
                tmpdir,
                "Output_of_taxonomy.csv"
            )
        
        results_path = tmpdir
        threads = 1
        accession_genus_dict = taxa_df.set_index("Genbank")["Genus"].to_dict()
        Figure = True
        verbose = True
        blastn_exe = "blastn"
        makeblastdb_exe = "makeblastdb"

        # Invoke the function under test
        merged_df, query_merged_df, closest_genome = classification_viridic(
            known_taxa_path,
            query,
            taxa_df,
            taxa_csv_output_path,
            results_path,
            threads,
            accession_genus_dict,
            Figure,
            verbose,
            blastn_exe,
            makeblastdb_exe
        )

        # Perform assertions
        assert isinstance(merged_df, pd.DataFrame)
        assert isinstance(query_merged_df, pd.DataFrame)
        assert isinstance(closest_genome, str)


    # Generates heatmaps and saves them if Figure is True
    def test_do_not_generates_heatmaps_if_figure_is_false(self, tmpdir):
        # Initialize input variables
        test_file_dir = os.path.dirname(__file__)
        
        # Initialize input variables
        known_taxa_path = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "known_taxa.fa"
                                 )
        query = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "query.fasta"
                                 )
        
        taxa_df = read_VMR(
            os.path.join(test_file_dir,
                        "dummy_genome",
                        "VMR.dummy.xlsx"
                    )
            )
        
        taxa_df = taxa_df.rename(columns={"Accession": "Genbank"})

        taxa_csv_output_path = os.path.join(
                tmpdir,
                "Output_of_taxonomy.csv"
            )
        
        results_path = tmpdir
        threads = 1
        accession_genus_dict = taxa_df.set_index("Genbank")["Genus"].to_dict()
        Figure = False
        verbose = True
        blastn_exe = "blastn"
        makeblastdb_exe = "makeblastdb"

        # Mock the necessary functions and objects

        # Invoke the function under test
        merged_df, query_merged_df, closest_genome = classification_viridic(
            known_taxa_path,
            query,
            taxa_df,
            taxa_csv_output_path,
            results_path,
            threads,
            accession_genus_dict,
            Figure,
            verbose,
            blastn_exe,
            makeblastdb_exe
        )

        # Perform assertions
        assert isinstance(merged_df, pd.DataFrame)
        assert isinstance(query_merged_df, pd.DataFrame)
        assert isinstance(closest_genome, str)

    def test_good_dataframe(self, tmpdir):
        # Initialize input variables
        test_file_dir = os.path.dirname(__file__)
        
        # Initialize input variables
        known_taxa_path = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "known_taxa.fa"
                                 )
        query = os.path.join(test_file_dir,
                                 "dummy_genome",
                                 "Results_per_genome",
                                 "MZ130489",
                                 "query.fasta"
                                 )
        
        taxa_df = read_VMR(
            os.path.join(test_file_dir,
                        "dummy_genome",
                        "VMR.dummy.xlsx"
                    )
            )
        
        taxa_df = taxa_df.rename(columns={"Accession": "Genbank"})

        taxa_csv_output_path = os.path.join(
                tmpdir,
                "Output_of_taxonomy.csv"
            )
        
        results_path = tmpdir
        threads = 1
        accession_genus_dict = taxa_df.set_index("Genbank")["Genus"].to_dict()
        Figure = False
        verbose = True
        blastn_exe = "blastn"
        makeblastdb_exe = "makeblastdb"

        # Mock the necessary functions and objects

        # Invoke the function under test
        merged_df, query_merged_df, closest_genome = classification_viridic(
            known_taxa_path,
            query,
            taxa_df,
            taxa_csv_output_path,
            results_path,
            threads,
            accession_genus_dict,
            Figure,
            verbose,
            blastn_exe,
            makeblastdb_exe
        )

        # Perform assertions
        expected_res = pd.read_table(
            os.path.join(test_file_dir,
                "dummy_genome",
                "Results_per_genome",
                "MZ130489",
                "Output_of_taxonomy.tsv"
            )
        )

        query_merged_df = query_merged_df.astype(str)
        expected_res = expected_res.astype(str)

        assert query_merged_df.equals(expected_res)


####################################################################################################