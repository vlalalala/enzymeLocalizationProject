import sys

def hello_function(name:str, exclamation_mark: bool):
    suffix = "!" if exclamation_mark else "."
    with open(f"{name}.txt", "a") as f:
        f.write(f"{name} was greeted{suffix}")

if __name__ == "__main__":
    print("Arguments received:", sys.argv)
    print(sys.argv[1])
    name = sys.argv[1]
    exclamation_mark = sys.argv[2].lower() == "true"
    hello_function(name, exclamation_mark)