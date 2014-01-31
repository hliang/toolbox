Collections of small and useful programs and scripts. Found on Internet or modified/created by me.

==

#### Create a new repository on the command line
```
touch README.md
git init
git add README.md
git commit -m "first commit"
git remote add origin https://github.com/hliang/toolbox.git
git push -u origin master
```

#### Push an existing repository from the command line
```
git remote add origin https://github.com/hliang/toolbox.git
git push -u origin master
```


#### After local repo files changed, pushes all the modified local objects to the remote repository and advances its branches.
`git push [alias] [branch]` `git push [remote-name] [branch-name]`.
```
git push origin master
```

you run `git push [alias] [branch]` to update a remote repository with the changes you've made locally. It will take what your [branch] looks like and push it to be [branch] on the remote, if possible.

If someone else has pushed since you last fetched and merged, the Git server will deny your push until you are up to date.
```
git fetch origin/master
git merge origin/master
```
`git fetch origin` -- Fetches all the objects from the remote repository that are not present in the local one.
`git merge newbranchversion` -- Merges one or more branches into your current branch and automatically creates a new commit if there are no conflicts.


#### If the remote origin repo changed, update your local repo.

Fetches the files from the remote repository and merges it with your local one. This command is equal to the `git fetch` and the `git merge` sequence.
```
git pull origin/upstream
```

