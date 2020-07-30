# How to contribute code using git
<small><i>Page last updated on: {{ git_revision_date }}</i></small>

## Part 1. One-time set up

### Identity:
```
$ git config --global user.name "<Your Name>"  
$ git config --global user.email <your email>
```

### Save credential for 6 hours:
```
$ git config --global credential.helper 'cache --timeout=21600'
```

### Save credential permanently:
```
$ git config --global credential.helper store
```

### Alias of showing git history:
Add this to ~/.gitconfig

```
[alias]
lg = log --graph --abbrev-commit --decorate --format=format:'%C(bold blue)%h%C(reset) - %C(bold green)(%ar)%C(reset) %C(white)%s%C(reset) %C(dim white)- %an%C(reset)%C(bold yellow)%d%C(reset)' --all
```

## Part 2. Start a new repository

![git picture](img/git.png)

### Clone a remote repo:
Clone a remote repo to local computer:
```
$ git clone <url>
```

### Convert a local directory:
Convert an existing local directory to local git repository:
```
$ git init
```

## Part 3. Single line development:

blockdiag {
default_fontsize = 20;

C1 [shape = circle];
C2 [shape = circle];
C3 [shape = circle];
C4 [shape = circle];

C1 -> C2;
C2 -> C3;
C3 -> C4;
}


### Add all files to track:
```
$ git add .
```

### Check status:
```
$ git status
```

### Commit:
```
$ git commit -a -m "<commit message>"
```

### Add a file to track:
```
$ git add <file name>
```

### Show history:
```
$ git lg
```
*[It requires you set up the “lg” alias in ~/.gitconfig](#alias-of-showing-git-history)*

## Part 4. Go to past commits

**Create a new branch and revert to a past commit:**
```
$ git checkout -b <new branch> <commit hash>
```

```<new branch>``` is a branch name you make.

```<commit hash>``` is the string as reported by git log

If you want to make this branch as the new master branch, do a swap as following:

1. make sure your are in the new branch
```
$ git checkout <new branch>
```

2. force master to merge with current branch and use current branch as favored:
```
$ git merge -s ours master
```

3. go to the master branch and reconcile again:
```
$ git checkout master
$ git merge <new branch>
```

4. after merge, delete the branch.
```
$ git branch -d <new branch>
```

## Part 5. Multi-line development
### Show branches:
```
$ git branch
```

###
Show remote branches:
```
$ git branch -r
```

### Switch between branches:
```
$ git checkout <branch name>
```

### Create and switch branch:
```
$ git checkout -b <branch name>
```

### Merge branches:
```
$ git merge <another branch>
```

This will merge a branch to current branch:

### Merge with current favored:
```
$ git merge -s ours <another branch>
```

### Delete branches
```
$ git branch -d <branch name>
```

### Force delete:
```
$ git branch -D <branch name>
```

## Part 6. Sync remote and local

### Check remote:
```
$ git remote -v
```
“origin” is the default name of your first remote.

### Add more remotes:
```
$ git remote add <remote> <url>
```

### Pull from remote repo:
```
$ git pull
```
Or
```
$ git pull <remote> <branch>
```

### Push to remote:
```
$ git push
```

### Push new local branch to remote:
```
$ git push -u <remote> <branch>
```

### Delete remote branch:
```
$ git push <remote> --delete <branch>
```

### Delete remote tracking branch:
```
$ git remote prune <remote>
```

## Part 7. Merge branches
### Common scenario of merge:
* Start a new feature:
```
$ git checkout -b new-feature
```

* Edit some files, then commit the change:
```
$ git commit -a -m "Start a feature"  
```

* Edit some files, then commit more changes:
```
$ git commit -a -m "Finish a feature"  
```

* Merge in the new-feature branch:
```
$ git checkout master 
$ git merge new-feature 
$ git branch -d new-feature
```

### Conflict in merge:
When conflicts occur, the conflicting files will have visual marks like:
```
<<<<<<< master
conflicting text in receiving branch
=======
conflicting text in merging branch
>>>>>>> branch
```

You need to edit text and remove <<<<<<, ======, >>>>>> lines.

Then run a commit:
```
$ git commit -a -m "<commit message>"
```