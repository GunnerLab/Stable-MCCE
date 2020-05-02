# How to edit MCCE wiki
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

This documentation is written in mkdocs, which uses markdown language and publishes on github pages at url https://github.com/GunnerLab/Stable-MCCE.

## Install of Mkdocs

You need pip, a python package manager, for installation. Generally, you are able to install pip through conda
```
conda install pip
```
or system package manager
```
sudo aptitude install python-pip
```

If you don't know pip, please go to https://pip.pypa.io/en/stable/ for more information.

In a terminal window, run these commands to install mkdocs, theme, and necessary plugins for documentation of this project.
```
pip install --upgrade pip
pip install bs4
pip install unicode
pip install mkdocs
pip install mkdocs-material
pip install pymdown-extensions
pip install markdown-blockdiag
pip install markdown-include
pip install mkdocs-git-revision-date-plugin
```

##  Writing and editing

All documentation files reside under folder wiki/

* mkdocs.yml : configuration file, also defines the menu
* docs/users/ : location of documentation files for users
* docs/developers/ : location of documentation files for developers
* site/ : site web pages published by command mkdocs gh-deploy

Under wiki/ directory, run
```
mkdocs serve
```

will bring up a web service, and point browser to http://localhost:8000 to enable viewing and debugging the documentation.

## Deploying

### Commit changes:
```
git commit -a -m "updated by jmao"
```

This commit will update the time stamp on the page as well.

### Push to your forked repository
First check if you set up the GunnerLab as upstream

For example:
```
(base) jmao@pc:~/projects/Stable-MCCE/wiki$ git remote -v
origin	https://github.com/newbooks/Stable-MCCE.git (fetch)
origin	https://github.com/newbooks/Stable-MCCE.git (push)
upstream	https://github.com/GunnerLab/Stable-MCCE.git (fetch)
upstream	https://github.com/GunnerLab/Stable-MCCE.git (push)
```

If you do not have upstream repo, push to the origin, then make a pull request:
```
git push
```

Then go to https://github.com/ under your forked repository. Make a pull request. Someone will publish the edit for you.

###  Publish directly
If you are a group member with publish permission, you can publish directly.

Make sure the upstream is setup properly:
```
git remote -v
```
If upstream is not pointing to https://github.com/GunnerLab/Stable-MCCE.git, run
```
git remote add upstream https://github.com/GunnerLab/Stable-MCCE.git
```

Under wiki folder

To deploy site, you can

 1. create a pull request so that the GunnerLab repo administrator can merge and deploy for you.

or

deploy directly by running gh-deploy under directory wiki/. You need GunnerLab repo member permission to do this.
If gh-deploy reports conflict in special branch gh-pages, use "--force" switch to overwrite:


mkdocs gh-deploy --force
The published site is on: https://gunnerlab.github.io/Develop-MCCE/.