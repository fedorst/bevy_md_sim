import subprocess
from pathlib import Path

# extensions to include
exts = {".json", ".md", ".ts", ".js", ".svelte", ".yml", ".env", ".py", ".rs", ".toml", ".py"}
extra_files = {"Dockerfile"}
out_path = Path("concatenated_dump.txt")

def get_git_included_files():
    result = subprocess.run(
        ["git", "ls-files", "--others", "--exclude-standard", "--cached"],
        capture_output=True, text=True
    )
    return [Path(line) for line in result.stdout.strip().splitlines()]

def is_included(path):
    return (
        (path.suffix in exts or path.name in extra_files)
        and path.name not in ["package-lock.json", "uv.lock", "next_up.ignore", "Cargo.lock", "dump_code_for_llm.py", "get_tree.sh"]
        and path.is_file()
    )

def get_file_tree():
    result = subprocess.run(
        ["./get_tree.sh"],
        capture_output=True, text=True
    )
    return result.stdout.strip()

def main():
    files = [f for f in get_git_included_files() if is_included(f)]
    files.sort()  # optional: deterministic order

    with out_path.open("w", encoding="utf-8") as out_file:
        out_file.write("/* -------- filetree -------- */\n")
        out_file.write(get_file_tree())
        out_file.write("\n\n")

        for f in files:
            out_file.write(f"/* ---- file: {f} ---- */\n")
            try:
                content = f.read_text(encoding="utf-8")
            except UnicodeDecodeError:
                content = f.read_text(encoding="latin1")  # fallback
            out_file.write(content)
            out_file.write("\n\n")

if __name__ == "__main__":
    main()
