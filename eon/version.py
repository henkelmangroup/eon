import subprocess
from os.path import abspath, dirname, join, exists

# this is old svn version
# def version():
#     try:
#         output = subprocess.check_output(["svnversion",str(path)], stderr=subprocess.STDOUT)
#     except subprocess.CalledProcessError as grepexc:                                                                                                   
#         return 'unknown'
#     return 'svn revision %s' % output.decode('ascii')

# version for git that will work even without pip install.
def version():
    path = dirname(abspath( __file__ ))
    main_dir = dirname(path)
    # first check if using release tarball (no .git folder)
    if not exists(join(main_dir,".git")):
        version_file = join(main_dir,".eon_release_version")
        if exists(version_file):
            with open(version_file) as f:
                return f.read().strip()
        return "unknown" # if no .git folder AND no version file, nothing can be done
    # check if exactly at a tag
    try:
        return subprocess.check_output(
            ["git", "describe", "--tags", "--exact-match"],
            cwd=path, stderr=subprocess.STDOUT
        ).decode("ascii").strip() # "v1.0.0" - exactly on a tag, done
    except subprocess.CalledProcessError:
        pass
    # if not on a tag, gather info needed for issue handling - last tag, branch, commit hash
    try:
        # if we cloned the repo and not at a tag, get the branch info
        branch = subprocess.check_output(["git","rev-parse","--abbrev-ref","HEAD"],stderr=subprocess.STDOUT,cwd=path)
        branch = branch.decode("ascii").strip()
        # now, get the specific commit for issue handling
        commit = subprocess.check_output(["git", "rev-parse","--short=8", "HEAD"],cwd=path,stderr=subprocess.STDOUT)
        commit = commit.decode("ascii").strip()
        # now, get the previous tag if it exists
        try:
            tag = subprocess.check_output(
                ["git", "describe", "--tags", "--abbrev=0"],
                cwd=path, stderr=subprocess.STDOUT
            ).decode("ascii").strip()
            return f"{tag} ({branch} #{commit})" # "v1.0.0 (main #abcdef12)"
        except subprocess.CalledProcessError:
            return f"{branch} #{commit}" # "main #abcdef12"
    except (subprocess.CalledProcessError, FileNotFoundError) as grepexc:
        return 'unknown'


if __name__ == '__main__':
    print(version())
