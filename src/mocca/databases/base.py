#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 15:53:30 2021

@author: haascp
"""

class BaseDatabase():
    def __init__(self):
        self.items = []
    
    def __iter__(self):
        """Yields all items inside the database.
           Allows statements such as "for component in component_database:" """
        for item in self.items:
            yield item
    
    def get_index(self, item):
        return self.items.index(item)
    
    def insert_item(self, item, unique=True):
        if unique:
            if item in self.items:
                idx = self.get_index(item)
                self.items[idx] = item
            else:
                self.items.append(item)
        else:
            self.items.append(item)
    
    def delete_all_items(self):
        """
        Clears all components out of database.
        """
        self.items = []
    
    def update_item(self, idx, item):
        self.items[idx] = item
    
    def get_item_by_index(self, idx):
        return self.items[idx]
    
    def get_next_obj_by_key_value(self, key, value):
        next(item for item in self.items if
                    item.key == value)
