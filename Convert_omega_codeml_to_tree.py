#!/usr/bin/env python3
import re

def convert_support_to_branchlength(tree_str):
    # Replace #number with :number, capping >1 at 1.0
    def replace_support(match):
        val = float(match.group(1))
        if val > 1:
            val = 1.0
        return f":{val}"
    return re.sub(r"#([\d\.]+)", replace_support, tree_str)

# Ask user to paste tree
print("Paste your Newick tree (end with Enter):")
tree_str = input().strip()

converted = convert_support_to_branchlength(tree_str)

print("\n✅ Converted tree:")
print(converted)
