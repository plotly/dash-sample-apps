Stuff I use with `plotly.py` and `dash` webapps.
Generally to be included as a submodule in a git repo.

To add to your repo:

```bash
git submodule add git@github.com:nicholas-esterer/plotly-common.git
```

Then after cloning your parent project, initialize the submodule (clone the repo
and checkout the commit registered in your parent project).

```bash
git submodule init
git submodule update
```
