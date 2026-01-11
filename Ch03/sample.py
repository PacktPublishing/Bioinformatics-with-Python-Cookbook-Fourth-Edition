#!/usr/bin/env python3

import os, sys  # Unused imports and multiple imports on one line


def example_function(a, b):  # Missing docstring
    if a > b:
        print("a is greater than b")  # Improper indentation
    else:
        print("b is greater or equal to a")  # Extra indentation


class ExampleClass:  # Missing class docstring
    def __init__(self, value):
        self.value = value  # Missing spaces around '=' operator
        self.data = []  # Unused attribute

    def add_data(self, item):  # Unused method argument 'item'
        pass

    def display(self):
        print("Value: ", self.value)  # Space before comma is bad style


# Unused variable and name not in snake_case
BADVariableName = 42

# Long line exceeding 80 characters
print(
    "This is a really, really, really, really, really, really, really long line of code."
)

example_function(10, 5)  # Function call with no meaningful context
