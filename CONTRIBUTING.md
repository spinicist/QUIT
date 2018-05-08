# Contributing to QUIT

Welcome to QUIT. Any help is very much appreciated.

These guidelines exist to make it easy to get involved and ensure that QUIT is welcoming to everyone.

## Asking Questions / Reporting Bugs

If you have general questions about how to use QUIT programs, you can either open an [issue](https://github.com/spinicist/QUIT/issues) or ask a question on (Neurostars)[https://neurostars.org].

If you think you have found a bug, then opening a Github issue is the preferred avenue. Please check whether an identical or similar issue already exists first. When opening a new issue, please give as much information as you can, including
- The version of QUIT you are using. Either give the version number if you downloaded the binaries, or a branch/git commit id if you compiled from source.
- The operating system you are running on.
- A description of the bug, including the full shell command you ran and any input files.

## Contributing Changes

If you want to edit the QUIT code yourself to fix a problem or implement a new feature, you are very welcome to! Please follow this model for submitting a Pull Request:
- (Fork)[https://help.github.com/articles/fork-a-repo/] the QUIT repo to your Github profile.
- Clone this repo, and checkout the `development` branch, then create a new branch.
- Make your changes on this own branch.
- When you are finished, make sure your copy of the `development` branch is up to date, and if necessary `rebase` your branch to the latest `development` commit.
- Push your branch to Github.
- Open a (Pull Request)[https://help.github.com/articles/fork-a-repo/].

To date, QUIT has one main developer (@spinicist). Hence the above is a suggested model for how changes can be contributed, and will be edited if and when a better model is suggested.

QUIT currently has a fairly loose coding style, in part because it uses several libraries with conflicting coding style. @spinicist prefers the following style at present, but will allow some flexibility (except as noted in the first point):
- Kernigan & Ritchie braces convention. Please do not use the GNU Coding Style braces convention as demonstrated in ITK.
- Classes and global functions should be Capitalised. Local variables should begin with a lower-case letter. Private/protected member variables should begin with `m_`, public member variables do not have to follow this to make the interface clearer.
- Camel case or `under_scores` are both acceptable as your whims dictate.
- Spaces, not tabs. If you find some tabs in one of the older files please let @spinicist know and we can decide what to do about them.

## Code of Conduct

- Anyone who participates in the development of QUIT is expected to show respect and courtesy to all other community members at all times.
- Harrasment in any form towards any members will not be tolerated.
- All communication should be appropraite for a professional audience including people of different backgrounds.
- Be kind to others. Do not insult or put down other contributers.

## Thanks

Above all, thank you for using or contributing to QUIT. It really is appreciated.
