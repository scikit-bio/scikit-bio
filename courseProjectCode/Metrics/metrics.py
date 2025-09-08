import os, ast
from pathlib import Path

def collect_metrics(directory):
    print("Maintainability Metrics")

    total_lines_of_code = 0
    total_comment_lines = []
    total_comment_density = []

    for py_file in sorted(Path(directory).rglob("*.py")):
        try:
            with open(py_file, "r", encoding="utf-8") as f:
                lines = f.readlines()
            
            total_lines = len(lines)
            total_lines_of_code += total_lines
            
            comment_lines = sum(1 for line in lines if line.strip().startswith("#"))
            comment_density = (comment_lines / total_lines * 100) if total_lines > 0 else 0

            total_comment_lines.append(comment_lines)
            total_comment_density.append(comment_density)
            
            print(f"{py_file.name}: {total_lines} LOC, {comment_density:.1f}% comments")

        except Exception as e:
            print(f"Error reading {py_file}: {e}")
    
    print("Testability Metrics")

    test_files = list(Path(directory).rglob("test_*.py")) + list(Path(directory).rglob("*_test.py"))
    total_test_functions = 0
    
    for test_file in sorted(test_files):
        try:
            with open(test_file, "r", encoding="utf-8") as f:
                tree = ast.parse(f.read())
                
            test_funcs = sum(1 for node in ast.walk(tree) if isinstance(node, ast.FunctionDef) and node.name.startswith("test_"))
            total_test_functions += test_funcs
            print(f"{test_file.name}: {test_funcs} test functions")
            
        except Exception as e:
            print(f"Error analyzing {test_file}: {e}")
    
    print(f"\nTotal test files: {len(test_files)}")
    print(f"Total test functions: {total_test_functions}")

    print("\n\nSummary Metrics")
    print(f"Total lines of code: {total_lines_of_code}")
    print(f"Total comment lines: {sum(total_comment_lines)}")
    print(f"Mean comment density: {sum(total_comment_density) / len(total_comment_density):.1f}%")
    print(f"Total test files: {len(test_files)}")
    print(f"Total test functions: {total_test_functions}")

if __name__ == "__main__":
    collect_metrics("skbio/")