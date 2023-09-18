#################################################################
# FILE : patient_diagnosis_BST.py
# WRITER : israel_nankencki , israelnan , 305702334
# EXERCISE : intro2cs2 ex11 2022
# DESCRIPTION: several classes for diagnose illness by its symptoms.
# STUDENTS I DISCUSSED THE EXERCISE WITH:none.
# WEB PAGES I USED: none
# NOTES:
#################################################################

##############################################################################
#                                   Imports                                  #
##############################################################################
import itertools
from typing import List, Dict, Union
from random import randint


PATH_LIST = List[bool]


class Node:
    """
    A class with binary nodes to build the diagnosis tree.
    """
    def __init__(self, data, positive_child=None, negative_child=None):
        self.data = data
        self.positive_child = positive_child
        self.negative_child = negative_child

    def node_minimiser(self, remove_empty=False) -> bool:
        """
        this is a recursive method that minimising diagnosis trees for Diagnose class.
        :param remove_empty: a bool flag for minimising a branch.
        :return: True if this branch need to be removed, False otherwise.
        """
        if not self.negative_child and not self.positive_child:
            return remove_empty
        if self.negative_child.node_minimiser():
            self.negative_child = self.negative_child.negative_child
        if self.positive_child.node_minimiser():
            self.positive_child = self.positive_child.positive_child
        if self.positive_child == self.negative_child:
            return True


class Record:
    """
    A class with record objects that have the name of any illness as a string, and list of strings for its symptoms.
    """
    def __init__(self, illness, symptoms):
        self.illness = illness
        self.symptoms = symptoms


def parse_data(filepath: str) -> List[Record]:
    """
    this function makes a list of Record object out of given text file.
    :param filepath: string with the text file path.
    :return: list of records objects.
    """
    with open(filepath) as data_file:
        records = []
        for line in data_file:
            words = line.strip().split()
            records.append(Record(words[0], words[1:]))
        return records


class Diagnoser:
    """
    A class for building and managing the diagnosis binary trees.
    """
    def __init__(self, root: Node) -> None:
        self.root = root

    def diagnose_helper(self, symptoms: List[str], node: Node) -> str:
        """
        this is a recursive helper method, which helps 'diagnose' to go over the diagnosis tree and diagnose.
        :param symptoms: list of reported symptoms as strings.
        :param node: the root of the diagnosis tree.
        :return: the name of the illness diagnosed as a string.
        """
        if not node.positive_child and not node.negative_child:
            return node.data
        if node.positive_child and node.data in symptoms:
            return self.diagnose_helper(symptoms, node.positive_child)
        if node.negative_child:
            return self.diagnose_helper(symptoms, node.negative_child)

    def diagnose(self, symptoms: List[str]) -> str:
        """
        this method diagnose the illness where it has been given a list of reported symptom.
        :param symptoms: list of reported symptoms as strings.
        :return: the name of the illness diagnosed as a string.
        """
        return self.diagnose_helper(symptoms, self.root)

    def calculate_success_rate(self, records: list[Record]) -> float:
        """
        this method calculates the success rate of diagnosis by checking it against a given list of Record objects.
        :param records: list of Record objects.
        :return: a float number for success rate of the diagnosis.
        """
        if not records:
            raise ValueError('records list has to contain some records.')
        correct_diagnose_score = 0
        for i in records:
            if self.diagnose(i.symptoms) == i.illness:
                correct_diagnose_score += 1
        return correct_diagnose_score / len(records)

    def _all_illnesses_helper(self, node: Node, all_illnesses_dict: Dict):
        """
        this is a recursive method that helps 'all_illnesses' to build a dictionary with illnesses name as key and its
         prevalence.
        :param node: the root of the diagnosis tree.
        :param all_illnesses_dict: dictionary with illnesses name as key and its prevalence.
        :return: none
        """
        if not node.positive_child and not node.negative_child:
            if node.data not in all_illnesses_dict:
                all_illnesses_dict[node.data] = 1
            else:
                all_illnesses_dict[node.data] += 1
        else:
            self._all_illnesses_helper(node.positive_child, all_illnesses_dict)
            self._all_illnesses_helper(node.negative_child, all_illnesses_dict)

    def all_illnesses(self) -> List[str]:
        """
        this method builds a list with all illness names contained in the diagnosis tree, sorted by its prevalence.
        :return: list with all illness names contained in the diagnosis tree, sorted by its prevalence.
        """
        all_illnesses_dict = {}
        self._all_illnesses_helper(self.root, all_illnesses_dict)
        all_illnesses_list = []
        for (ill, key) in sorted(all_illnesses_dict.items(), key=lambda x: x[1]):
            if ill:
                all_illnesses_list.append(ill)
        return all_illnesses_list[::-1]

    def _paths_to_illness_helper(self, node: Node, illness: Union[str, None]) -> List[PATH_LIST]:
        """
        this is a recursive method that helps 'paths_to_illness' to go over the diagnosis tree and build a list with
         paths to a given illness.
        :param node: the root of the diagnosis tree.
        :param illness: name of the given illness.
        :return: list with lists of paths to the given illness, contain True if it went through positive child or False
         otherwise .
        """
        if not node.positive_child and not node.negative_child:
            if node.data == illness:
                return [[]]
            return []
        all_paths = []
        for positive_path in self._paths_to_illness_helper(node.positive_child, illness):
            positive_path.insert(0, True)
            all_paths.append(positive_path)
        for negative_path in self._paths_to_illness_helper(node.negative_child, illness):
            negative_path.insert(0, False)
            all_paths.append(negative_path)
        return all_paths

    def paths_to_illness(self, illness: Union[str, None]) -> List[PATH_LIST]:
        """
        this method goes over the diagnosis tree and builds a list with paths to a given illness.
        :param illness: name of the given illness.
        :return: list with lists of paths to the given illness, contain True if it went through positive child or False
         otherwise .
        """
        return self._paths_to_illness_helper(self.root, illness)

    def minimise(self, remove_empty=False):
        """
        this method minimizes the diagnosis tree to a more efficient tree.
        :param remove_empty: a bool flag for minimising a branch.
        :return: the more efficient diagnosis tree, without the branches that has no effect on the final diagnosis or
         questions in the way to it.
        """
        Node.node_minimiser(self.root, remove_empty)
        return self


def sort_frequency_dict_to_list(all_illnesses_dict: Dict[str, int]) -> List[str]:
    """
    this function helps 'get_match_or_random' to create a list of all illnesses in the tree, sorted by their frequency.
    :param all_illnesses_dict: dictionary with all illness names in the tree as keys, and its frequency as value.
    :return: list of all illnesses in the tree, sorted by their frequency.
    """
    return sorted(all_illnesses_dict, key=all_illnesses_dict.get, reverse=True)


def is_it_the_same_symptom(record: Record, negative_symp: List[str]) -> bool:
    """
    this function helps 'get_match_or_random' to whether one of the negative path of symptoms in the tree is in the
     record list of symptoms.
    :param record: Record object to compare with.
    :param negative_symp: list of symptoms as strings for the negative path through the tree.
    :return: True if the at least one of the symptoms is in the record symptoms list, False otherwise.
    """
    for symp in negative_symp:
        if symp in record.symptoms:
            return True
    return False


def is_it_match(record: Record, positive_symptoms: List[str]) -> bool:
    """
    this function helps 'get_match_or_random' to whether the path of symptoms in the tree is correlated with its record.
    :param record: Record object to compare with.
    :param positive_symptoms: list of symptoms as strings for the path through the tree to an illness.
    :return: True if the path to an illness in the tree is fully correlated with its record, False otherwise.
    """
    for symp in positive_symptoms:
        if symp not in record.symptoms:
            return False
    return True


def get_match_or_random(records: List[Record], negative_symp: List[str], positive_symp: List[str]) -> str:
    """
    this function helps '_build_tree_helper' to crate the final leafs of the diagnosis tree.
    :param records: list of Record objects.
    :param negative_symp: list of symptoms that the way to it is through negative child in the diagnosis tree.
    :param positive_symp: list of symptoms that the way to it is through positive child in the diagnosis tree.
    :return: the illness name for the given leaf as string.
    """
    all_illnesses_dict = {}
    for rec in records:
        if is_it_match(rec, positive_symp) and not is_it_the_same_symptom(rec, negative_symp):
            if rec.illness not in all_illnesses_dict:
                all_illnesses_dict[rec.illness] = 1
            else:
                all_illnesses_dict[rec.illness] += 1
    if all_illnesses_dict:
        return sort_frequency_dict_to_list(all_illnesses_dict)[0]
    return records[randint(0, len(records) - 1)].illness


def _build_tree_helper(records: List[Record], symptoms: List[str], negative_symp: List[str],
                       positive_symp: List[str]) -> Diagnoser:
    """
    this is a recursive function that helps 'build_tree' to build a Diagnoser object with diagnosis tree out of a
     given lists of Record objects and symptoms.
    :param records: list of Record objects.
    :param symptoms: list of symptoms as strings.
    :param negative_symp: list of symptoms that the way to it is through negative child in the diagnosis tree.
    :param positive_symp: list of symptoms that the way to it is through positive child in the diagnosis tree.
    :return: Diagnoser object with diagnosis tree in it.
    """
    if not symptoms:
        return Diagnoser(Node(get_match_or_random(records, positive_symp, negative_symp)))
    positive_child = _build_tree_helper(records, symptoms[1:], negative_symp, positive_symp + [symptoms[0]])
    negative_child = _build_tree_helper(records, symptoms[1:], negative_symp + [symptoms[0]], positive_symp)
    return Diagnoser(Node(symptoms[0], positive_child, negative_child))


def is_it_only_records_and_strings(records: List[Record], symptoms: List[str]) -> bool:
    """
    this function helps 'build_tree' and 'optimal_tree' to check whether all objects within their given lists of
     symptoms and records are legal, and raise an exception if needed.
    :param records: list of Record objects.
    :param symptoms: list of symptoms as strings.
    :return: True if both lists are legal, no return otherwise.
    """
    for rec in records:
        if type(rec) != Record:
            raise TypeError("records list must contain Record type objects only")
    for symp in symptoms:
        if type(symp) != str:
            raise TypeError("symptoms list must contain string objects only")
    return True


def build_tree(records: List[Record], symptoms: List[str]) -> Diagnoser:
    """
    this function builds Diagnoser object with diagnosis tree out of a given lists of Record objects and symptoms.
    :param records: list of Record objects.
    :param symptoms: list of symptoms as strings.
    :return: Diagnoser object with diagnosis tree in it.
    """
    if not records:
        return Diagnoser(Node(None))
    if is_it_only_records_and_strings(records, symptoms):
        return _build_tree_helper(records, symptoms, [], [])


def optimal_tree(records: List[Record], symptoms: List[str], depth: int) -> Diagnoser:
    """
    this function builds an optimal diagnosis tree with the higher success rate among all combinations of symptoms'
     subgroups.
    :param records: list of Record objects.
    :param symptoms: list of symptoms as strings.
    :param depth: int for the maximum size of the subgroups.
    :return: Diagnoser object with the most optimal diagnosis tree in it.
    """
    if is_it_only_records_and_strings(records, symptoms):
        success_rate = 0
        optimal_diagnoser_tree = Diagnoser(Node(None))
        for symp_comb in itertools.combinations(symptoms, depth):
            cur_root = build_tree(records, symp_comb)
            cur_root_success_rate = cur_root.calculate_success_rate(records)
            if cur_root_success_rate > success_rate:
                success_rate, optimal_diagnoser_tree = cur_root_success_rate, cur_root
        return optimal_diagnoser_tree


if __name__ == "__main__":

    # Manually build a simple tree.
    #                cough
    #          Yes /       \ No
    #        fever           healthy
    #   Yes /     \ No
    # covid-19   cold

    flu_leaf = Node("covid-19", None, None)
    cold_leaf = Node("cold", None, None)
    inner_vertex = Node("fever", flu_leaf, cold_leaf)
    healthy_leaf = Node("healthy", None, None)
    root = Node("cough", inner_vertex, healthy_leaf)
    diagnoser = Diagnoser(root)

    # Simple test
    diagnosis = diagnoser.diagnose(["cough"])
    if diagnosis == "cold":
        print("Test passed")
    else:
        print("Test failed. Should have printed cold, printed: ", diagnosis)
