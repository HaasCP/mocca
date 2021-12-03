#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  2 15:53:30 2021

@author: haascp
"""

class BaseDatabase():
    """
    Base class for all databases with the shared primary key compound_id.
    """
    def __init__(self):
        self.items = []
    
    def __iter__(self):
        """Yields all items inside the database.
           Allows statements such as "for component in component_database:" """
        for item in self.items:
            yield item

    def __contains__(self, compound_id):
        """Allows statements such as "'Component 1' in component_database"
        to see if that Component is inside the database"""
        return compound_id in [item.compound_id for item in self.items]

    def __getitem__(self, compound_id):
        """Allows statements such as "component_database['Component 1']" to access
           that Component"""
        if compound_id not in self:
            return next(item for item in self.items if
                        item.compound_id == compound_id)
        else:
            raise AttributeError("{} not found in Database!"
                                 "".format(compound_id))

    def update_item(self, old_item, new_item):
        self.items = [new_item if item == old_item else item for item in self.items]

    def insert_item(self, item, unique=True):
        if unique:
            if item in self.items:
                raise NotImplementedError
            else:
                self.items.append(item)
        else:
            self.items.append(item)

    def delete_all_items(self):
        """
        Clears all components out of database.
        """
        self.items = []
