"""U25ADG02 - Data Structures and Algorithms | Group 7 Project
Biomedical Applications of Expression Trees, Graph Traversal and Algorithm Analysis
Modules: clinical expression trees, BFS contact tracing, recursive/iterative growth.
All contact-network data is synthetic. This is an educational demonstration, not a clinical tool.
Run with: python G07_project.py
Optional graph: pip install matplotlib
"""
import math
import time
import sys
sys.setrecursionlimit(5000)

# --------------------------- Stack and queue ---------------------------
class Stack:
    class Node:
        def __init__(self, value, nxt=None):
            self.value, self.next = value, nxt
    def __init__(self): self.top = None
    def push(self, value): self.top = self.Node(value, self.top)
    def pop(self):
        if self.top is None: raise IndexError("stack empty")
        value, self.top = self.top.value, self.top.next
        return value
    def peek(self): return None if self.top is None else self.top.value
    def is_empty(self): return self.top is None

class Queue:
    class Node:
        def __init__(self, value): self.value, self.next = value, None
    def __init__(self): self.head = self.tail = None
    def enqueue(self, value):
        node = self.Node(value)
        if self.tail is None: self.head = self.tail = node
        else: self.tail.next = node; self.tail = node
    def dequeue(self):
        if self.head is None: raise IndexError("queue empty")
        value = self.head.value
        self.head = self.head.next
        if self.head is None: self.tail = None
        return value
    def is_empty(self): return self.head is None

# --------------------------- Module 1: expression trees ---------------------------
PREC = {'+': 1, '-': 1, '*': 2, '/': 2, '^': 3}
RIGHT_ASSOC = {'^'}
FORMULAS = {
    '1': ("BMI = weight(kg) / height(m)^2", "w / h ^ 2"),
    '2': ("BSA (Mosteller) = sqrt(height(cm)*weight(kg)/3600)", "( h * w / 3600 ) ^ 0.5"),
    '3': ("MAP = (SBP + 2*DBP) / 3", "( s + 2 * d ) / 3"),
}

class TNode:
    def __init__(self, value, left=None, right=None):
        self.value, self.left, self.right = value, left, right

def tokenize(expression):
    tokens, i = [], 0
    while i < len(expression):
        c = expression[i]
        if c.isspace(): i += 1
        elif c.isdigit() or c == '.':
            j = i
            while j < len(expression) and (expression[j].isdigit() or expression[j] == '.'): j += 1
            token = expression[i:j]
            try: float(token)
            except ValueError: raise ValueError("Invalid number: " + token)
            tokens.append(token); i = j
        elif c.isalpha() or c == '_':
            j = i
            while j < len(expression) and (expression[j].isalnum() or expression[j] == '_'): j += 1
            tokens.append(expression[i:j]); i = j
        elif c in PREC or c in '()': tokens.append(c); i += 1
        else: raise ValueError("Invalid character: %s" % c)
    if not tokens: raise ValueError("Empty expression")
    return tokens

def infix_to_postfix(tokens):
    output, operators = [], Stack()
    for token in tokens:
        if token in PREC:
            while (not operators.is_empty() and operators.peek() in PREC and
                   (PREC[operators.peek()] > PREC[token] or
                    (PREC[operators.peek()] == PREC[token] and token not in RIGHT_ASSOC))):
                output.append(operators.pop())
            operators.push(token)
        elif token == '(': operators.push(token)
        elif token == ')':
            while not operators.is_empty() and operators.peek() != '(': output.append(operators.pop())
            if operators.is_empty(): raise ValueError("Mismatched brackets")
            operators.pop()
        else: output.append(token)
    while not operators.is_empty():
        token = operators.pop()
        if token == '(': raise ValueError("Mismatched brackets")
        output.append(token)
    return output

def build_tree(postfix):
    stack = Stack()
    for token in postfix:
        if token in PREC:
            if stack.is_empty(): raise ValueError("Malformed expression")
            right = stack.pop()
            if stack.is_empty(): raise ValueError("Malformed expression")
            left = stack.pop()
            stack.push(TNode(token, left, right))
        else: stack.push(TNode(token))
    if stack.is_empty(): raise ValueError("Empty expression")
    root = stack.pop()
    if not stack.is_empty(): raise ValueError("Malformed expression")
    return root

def inorder(node):
    if node.left is None and node.right is None: return node.value
    return '(' + inorder(node.left) + ' ' + node.value + ' ' + inorder(node.right) + ')'
def preorder(node):
    if node is None: return ''
    return (node.value + ' ' + preorder(node.left) + ' ' + preorder(node.right)).strip()
def postorder(node):
    if node is None: return ''
    return (postorder(node.left) + ' ' + postorder(node.right) + ' ' + node.value).strip()
def evaluate(node, values):
    if node.left is None and node.right is None:
        try: return float(node.value)
        except ValueError:
            if node.value not in values: raise KeyError("No value for variable '%s'" % node.value)
            return values[node.value]
    a, b = evaluate(node.left, values), evaluate(node.right, values)
    if node.value == '+': return a + b
    if node.value == '-': return a - b
    if node.value == '*': return a * b
    if node.value == '/':
        if b == 0: raise ZeroDivisionError("Division by zero")
        return a / b
    if node.value == '^':
        if a < 0 and not float(b).is_integer(): raise ValueError("Negative base with fractional power")
        result = a ** b
        if isinstance(result, complex): raise ValueError("Expression produced a complex result")
        return result
def tree_height(node): return 0 if node is None else 1 + max(tree_height(node.left), tree_height(node.right))
def print_tree(node, level=0):
    if node is None: return
    print_tree(node.right, level + 1)
    print('       ' * level + node.value)
    print_tree(node.left, level + 1)

def read_float(prompt, low=None, high=None):
    while True:
        try:
            value = float(input(prompt))
            if not math.isfinite(value): print("Enter a finite number."); continue
            if low is not None and value < low or high is not None and value > high:
                print("Value must be in the permitted range."); continue
            return value
        except ValueError: print("Invalid number; try again.")

def module1():
    while True:
        print("\n--- MODULE 1: Expression Tree for Clinical Formulas ---")
        print("1. BMI   2. BSA   3. MAP   4. Custom expression   0. Back")
        choice = input("Choice: ").strip()
        if choice == '0': return
        if choice in FORMULAS: title, expression = FORMULAS[choice]; print("Formula:", title)
        elif choice == '4': expression = input("Enter infix expression (+ - * / ^ and brackets; no unary minus): "); title = "Custom expression"
        else: print("Invalid choice."); continue
        try:
            tokens = tokenize(expression); postfix = infix_to_postfix(tokens); root = build_tree(postfix)
            print("Postfix :", ' '.join(postfix)); print("Inorder :", inorder(root))
            print("Preorder:", preorder(root)); print("Postorder:", postorder(root))
            print("Tree height:", tree_height(root)); print("Tree (rotated, root at left):"); print_tree(root)
            variables = []
            for token in tokens:
                if token not in PREC and token not in ('(', ')'):
                    try: float(token)
                    except ValueError:
                        if token not in variables: variables.append(token)
            env = {}
            for var in variables:
                lower = 0 if choice != '4' else None
                env[var] = read_float("Enter value for %s: " % var, lower)
            print("Result = %.4f" % evaluate(root, env))
        except (ValueError, ZeroDivisionError, KeyError, OverflowError) as error: print("Error:", error)

# --------------------------- Module 2: BFS contact tracing ---------------------------
class ContactGraph:
    def __init__(self): self.names, self.adj = [], []
    def index_of(self, name):
        for i, item in enumerate(self.names):
            if item == name: return i
        return -1
    def add_person(self, name):
        if not name or self.index_of(name) != -1: return False
        self.names.append(name); self.adj.append([]); return True
    def add_contact(self, first, second):
        i, j = self.index_of(first), self.index_of(second)
        if i == -1 or j == -1 or i == j or j in self.adj[i]: return False
        self.adj[i].append(j); self.adj[j].append(i); return True
    def bfs(self, source):
        levels, parents = [-1] * len(self.names), [-1] * len(self.names)
        queue = Queue(); levels[source] = 0; queue.enqueue(source)
        while not queue.is_empty():
            current = queue.dequeue()
            for neighbor in self.adj[current]:
                if levels[neighbor] == -1:
                    levels[neighbor] = levels[current] + 1; parents[neighbor] = current; queue.enqueue(neighbor)
        return levels, parents
    def chain(self, first, second):
        i, j = self.index_of(first), self.index_of(second)
        if i == -1 or j == -1: return None
        levels, parents = self.bfs(i)
        if levels[j] == -1: return None
        path, current = [], j
        while current != -1: path.append(self.names[current]); current = parents[current]
        return path[::-1]

def sample_graph():
    graph = ContactGraph()
    for number in range(1, 25): graph.add_person("P%02d" % number)
    edges = [(1,2),(1,3),(1,4),(2,5),(2,6),(3,6),(3,7),(4,8),(5,9),(6,10),(7,10),(7,11),(8,12),(9,13),(10,14),(11,15),(12,16),(13,17),(14,18),(15,19),(16,20),(17,21),(18,21),(22,23),(23,24)]
    for a, b in edges: graph.add_contact("P%02d" % a, "P%02d" % b)
    return graph

def module2(graph):
    while True:
        print("\n--- MODULE 2: BFS Contact Tracing ---")
        print("1. Show contact list  2. Trace from index patient  3. Shortest chain")
        print("4. Add person  5. Add contact  0. Back")
        choice = input("Choice: ").strip()
        if choice == '0': return
        elif choice == '1':
            for i, name in enumerate(graph.names):
                print("%s -> %s" % (name, ', '.join(graph.names[j] for j in graph.adj[i]) or '(none)'))
        elif choice == '2':
            source = graph.index_of(input("Index patient ID (e.g. P01): ").strip().upper())
            if source == -1: print("Unknown person."); continue
            levels, _ = graph.bfs(source)
            for level in range(max(levels) + 1):
                group = [graph.names[i] for i, value in enumerate(levels) if value == level]
                print("Level %d: %s" % (level, ', '.join(group)))
            missing = [graph.names[i] for i, value in enumerate(levels) if value == -1]
            print("Not connected:", ', '.join(missing) if missing else 'none')
        elif choice == '3':
            first = input("From: ").strip().upper(); second = input("To: ").strip().upper()
            if graph.index_of(first) == -1 or graph.index_of(second) == -1: print("Unknown person."); continue
            path = graph.chain(first, second)
            print("No chain: not connected." if path is None else "Chain (%d contacts): %s" % (len(path)-1, ' -> '.join(path)))
        elif choice == '4':
            name = input("New ID: ").strip().upper()
            print("Added." if graph.add_person(name) else "Invalid or duplicate ID.")
        elif choice == '5':
            first = input("Person A: ").strip().upper(); second = input("Person B: ").strip().upper()
            print("Contact added." if graph.add_contact(first, second) else "Invalid, same person, or contact already exists.")
        else: print("Invalid choice.")

# --------------------------- Module 3: bacterial growth ---------------------------
def growth_recursive(n, n0, rate, counter):
    counter[1] += 1
    if n == 0: return n0
    counter[0] += 1
    return rate * growth_recursive(n - 1, n0, rate, counter)
def growth_iterative(n, n0, rate, counter):
    value = n0
    for _ in range(n): value *= rate; counter[0] += 1
    return value
def growth_fast(n, n0, rate):
    result, base, exponent, multiplications = 1, rate, n, 0
    while exponent > 0:
        if exponent & 1: result *= base; multiplications += 1
        exponent >>= 1
        if exponent: base *= base; multiplications += 1
    return n0 * result, multiplications

def module3():
    while True:
        print("\n--- MODULE 3: Bacterial Growth ---")
        print("1. Compute colony count  2. Compare operation counts and timings  3. Recurrence analysis  0. Back")
        choice = input("Choice: ").strip()
        if choice == '0': return
        if choice == '1':
            n0 = read_float("Initial colony count N0 (1 to 1e9): ", 1, 1e9)
            rate = read_float("Growth factor r (1 to 10): ", 1, 10)
            n = int(read_float("Time steps n (0 to 1000): ", 0, 1000))
            rec, itr = [0,0], [0,0]
            try:
                a = growth_recursive(n, n0, rate, rec); b = growth_iterative(n, n0, rate, itr)
                print("Recursive: %.6g (multiplications=%d, calls=%d)" % (a, rec[0], rec[1]))
                print("Iterative: %.6g (multiplications=%d)" % (b, itr[0]))
            except (RecursionError, OverflowError) as error: print("Calculation too large:", error)
        elif choice == '2':
            sizes = [10, 50, 100, 200, 400, 800]
            print("n  recursive_mult  recursive_calls  iterative_mult  recursive_us  iterative_us")
            rec_counts, iter_counts = [], []
            for n in sizes:
                c1 = [0,0]; c2 = [0,0]
                start = time.perf_counter()
                for _ in range(100): c1 = [0,0]; growth_recursive(n, 1, 2, c1)
                rec_us = (time.perf_counter()-start) * 1e4
                start = time.perf_counter()
                for _ in range(100): c2 = [0,0]; growth_iterative(n, 1, 2, c2)
                iter_us = (time.perf_counter()-start) * 1e4
                rec_counts.append(c1[0]); iter_counts.append(c2[0])
                print("%-3d %-15d %-16d %-15d %-13.2f %.2f" % (n,c1[0],c1[1],c2[0],rec_us,iter_us))
            try:
                import matplotlib
                matplotlib.use('Agg')
                import matplotlib.pyplot as plt
                plt.plot(sizes, rec_counts, 'o-', label='Recursive multiplications')
                plt.plot(sizes, iter_counts, 's--', label='Iterative multiplications')
                plt.xlabel('Input size n (time steps)'); plt.ylabel('Operation count')
                plt.title('Bacterial growth: operations vs n'); plt.legend(); plt.grid(True)
                plt.savefig('growth_operations.png', dpi=150); plt.close()
                print("Graph saved as growth_operations.png")
            except ImportError: print("Optional graph skipped. Install with: pip install matplotlib")
        elif choice == '3':
            print("N(t) = r*N(t-1), N(0)=N0\nT(n)=T(n-1)+c, T(0)=d\nBackward substitution: T(n)=T(0)+n*c=c*n+d\nRecursive: Theta(n) time, Theta(n) stack space\nIterative: Theta(n) time, Theta(1) extra space\nExponentiation by squaring: Theta(log n) multiplications.")
        else: print("Invalid choice.")

# --------------------------- Tests ---------------------------
def run_tests():
    passed = total = 0
    def check(name, expected, actual):
        nonlocal passed, total
        total += 1; ok = expected == actual
        if ok: passed += 1
        print("%-4s %-40s expected=%s actual=%s" % ("PASS" if ok else "FAIL", name, expected, actual))
    print("\n=== Module 1 tests ===")
    root = build_tree(infix_to_postfix(tokenize("w / h ^ 2")))
    check("BMI w=70, h=1.75", "22.8571", "%.4f" % evaluate(root, {'w':70,'h':1.75}))
    expr = "( s + 2 * d ) / 3"; root = build_tree(infix_to_postfix(tokenize(expr)))
    check("MAP postfix", "s 2 d * + 3 /", ' '.join(infix_to_postfix(tokenize(expr))))
    check("MAP value", "93.3333", "%.4f" % evaluate(root, {'s':120,'d':80}))
    try: evaluate(build_tree(infix_to_postfix(tokenize("a/b"))), {'a':1,'b':0}); result = "no error"
    except ZeroDivisionError: result = "ZeroDivisionError"
    check("Division by zero", "ZeroDivisionError", result)
    try: infix_to_postfix(tokenize("(1+2")); result = "no error"
    except ValueError: result = "ValueError"
    check("Mismatched brackets", "ValueError", result)
    print("\n=== Module 2 tests ===")
    graph = sample_graph(); levels, _ = graph.bfs(graph.index_of("P01"))
    check("P06 level from P01", 2, levels[graph.index_of("P06")])
    check("Shortest chain P01 to P10", "P01 -> P02 -> P06 -> P10", ' -> '.join(graph.chain("P01","P10")))
    check("P22 unreachable", True, graph.chain("P01","P22") is None)
    check("Duplicate contact rejected", False, graph.add_contact("P01","P02"))
    print("\n=== Module 3 tests ===")
    counter = [0,0]; check("Recursive n=10", 1024, growth_recursive(10,1,2,counter))
    counter = [0,0]; actual = (growth_iterative(0,5,2,counter), counter[0]); check("Iterative n=0", (5,0), actual)
    c1, c2 = [0,0], [0,0]; a = growth_recursive(50,3,2,c1); b = growth_iterative(50,3,2,c2)
    check("Recursive equals iterative n=50", True, a == b and c1[0] == c2[0] == 50)
    check("Fast method n=30", 2**30, growth_fast(30,1,2)[0])
    print("\nPassed %d of %d tests" % (passed,total))

# --------------------------- Main menu ---------------------------
def main():
    graph = sample_graph()
    while True:
        print("\n========== Group 7 - Biomedical DSA Project ==========")
        print("1. Module 1 - Expression tree (clinical formulas)")
        print("2. Module 2 - BFS contact tracing")
        print("3. Module 3 - Bacterial growth analysis")
        print("4. Run all built-in tests")
        print("0. Exit")
        choice = input("Choice: ").strip()
        if choice == '1': module1()
        elif choice == '2': module2(graph)
        elif choice == '3': module3()
        elif choice == '4': run_tests()
        elif choice == '0': print("Goodbye."); break
        else: print("Invalid choice. Enter 1, 2, 3, 4 or 0.")

if __name__ == '__main__': main()
