from intervaltree import Interval, IntervalTree
from typing import Iterator, List, Callable
import copy
import os

from pyfaidx import Fasta

from ..classes.treader import TReader
from ..classes.transcript import Transcript, Object
from ..utils.common import *

class SpliceGraph:
    def __init__(self):
        self.nodes = []
        self.edges = []  # Will now store tuples of (target_node, edge_ids)

    def clear(self):
        self.nodes = []
        self.edges = []

    def add_node(self, node):
        self.nodes.append(node)
        self.edges.append([])
        return len(self.nodes) - 1

    def add_edge(self, from_node, to_node):
        # Find intersection of IDs between connected nodes
        edge_ids = set(self.nodes[from_node][2]) & set(self.nodes[to_node][2])
        self.edges[from_node].append((to_node, edge_ids))

    def add_from_chains(self, chains):
        # assumes intervals are sorted by start position
        id_map = {y: None for x in chains for y in x[2]}  # for each ID in the intervals, maps to the last inserted node
        for interval in chains:
            node_id = self.add_node(interval)
            for id in interval[2]:
                if id_map[id] is not None:
                    self.add_edge(id_map[id], node_id)
                id_map[id] = node_id

    def find_chain_path(self, target_chain, start_node_idx=0):
        """
        Check if target chain exists as a path in the graph by consuming pieces of the target.
        Can start at a specific node index to avoid redundant searches.
        Returns intersection of IDs collected along the valid path and the index of the first matching node.
        """
        if not target_chain or not self.nodes:
            return set()

        def explore_from_nodes(edges_to_try, remaining_target, current_ids):
            """Helper function to try path finding from multiple nodes"""
            # base case - if no remaining target, return current IDs
            if len(remaining_target) == 0:
                return current_ids
            for next_node, edge_ids in edges_to_try:
                # Update IDs to only include those that survive the edge transition
                next_ids = current_ids & edge_ids if current_ids is not None else edge_ids
                if next_ids:  # Only continue if we have surviving IDs
                    result = find_path_from_node(next_node, remaining_target, next_ids)
                    if result is not None:
                        return result
                    else: # dead end - remove this edge from the list
                        current_ids = current_ids - edge_ids
            return  None

        def find_path_from_node(node_idx, remaining_target, current_ids):
            """
            Recursive DFS that consumes pieces of the target interval.
            Returns None if path is invalid, set of IDs if valid path found.
            """
            # Base case - successfully consumed entire target
            if not remaining_target:
                return current_ids

            node = self.nodes[node_idx]
            node_start, node_end = node[0], node[1]
            target_start, target_end = remaining_target[0]
            
            # If node cannot overlap with target start, fail
            if node_end < target_start or node_start > target_end:
                return None
                
            # Update ID intersection
            new_ids = current_ids & set(node[2]) if current_ids is not None else set(node[2])
            if not new_ids:
                return None

            # If node fully contains target interval
            if node_start <= target_start and node_end >= target_end:
                # Continue with children nodes for the next target interval
                return explore_from_nodes(self.edges[node_idx], remaining_target[1:], new_ids)
            
            # If node covers start of target but not end
            if node_start <= target_start and node_end < target_end:
                new_target = [(node_end + 1, target_end)] + remaining_target[1:]
                # Try all possible next nodes
                return explore_from_nodes(self.edges[node_idx], new_target, new_ids)
            
            # If node starts after target start
            if node_start > target_start:
                # Need to find a path that covers the gap before this node
                new_target = [(target_start, node_start - 1)] + [(node_start, target_end)] + remaining_target[1:]
                # Try continuing from current node and its siblings
                return find_path_from_node(node_idx, new_target, new_ids)
                    
            return None

        # # Convert target chain to list of tuples for easier manipulation
        # target = [(interval[0], interval[1]) for interval in target_chain]
        
        # Try starting from each possible node
        first_start_node_idx = None
        for start_idx, node in enumerate(self.nodes[start_node_idx:], start_node_idx):
            if node[1] >= target_chain[0][0]:  # Node can cover start of target
                first_start_node_idx = start_idx if first_start_node_idx is None else first_start_node_idx # store the first node that can cover the target
                result = find_path_from_node(start_idx, target_chain, None)
                if result is not None:
                    return result, first_start_node_idx
            if node[0] > target_chain[0][1]: # can safely terminate search if we've passed the target
                break

        return set(), first_start_node_idx