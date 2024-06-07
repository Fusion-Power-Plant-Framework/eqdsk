# Contributing to ``eqdsk``

When contributing to this repository, please first discuss the change you wish to make
via issue, email, or any other method with the owners of this repository before making a
change.

## Code of Conduct

This project and everyone participating in it is governed by the  Contributor Covenant
Code of Conduct. By participating, you are expected to uphold this code. Please report
unacceptable behavior to [oliver.funk@ukaea.uk](mailto:olive.funk@ukaea.uk) and/or [james.cook1@ukaea.uk](mailto:james.cook1@ukaea.uk).


## Submitting an issue

In order to help us address your issue effectively, we ask that you follow this
procedure when creating an issue:

* Check that the issue is not already described elsewhere in [Issues
  ](https://github.com/Fusion-Power-Plant-Framework/eqdsk/issues)
* Write a fairly complete description of the bug/feature/problem in the issue, using
  the provided templates and adding the relevant tags
* If the issue is linked to a [Project](https://github.com/Fusion-Power-Plant-Framework/eqdsk/projects), please tag it

## Submitting a bug report

``Eqdsk`` is software in development and is therefore likely to contain bugs. If you
discover bugs, please follow this procedure:

* Raise an issue using the bug template with a `bug` flag
* Include a way to reproduce the bug in the issue, let us know the expected result and
  what actually happens

## Submitting a pull request

Please discuss any feature ideas you have with the developers before submitting them, as
you may not be aware of parallel development work taking place, or implementation
decisions / subtleties which are relevant. The ideal workflow for submitting a pull
request is as follows:

* Discuss the feature with the developers first
* Submit an issue documenting the intent of the feature you wish to develop
* Fork our repository (see the [GitHub documentation on forking](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks))
* Make a branch of the `develop` branch which bears a similar name as your issue (e.g.
  `new_feature`)
* Develop your feature(s) in your `new_feature` branch
* Discuss any problems that arise in the associated issue, perhaps documenting your
  progress
* Finally, as the author of the `new_feature` branch, you must submit a pull request
  onto the `develop` branch
  * Link the relevant issue(s) and project(s)
  * Add an assignee, who will be responsible for reviewing and merging the PR. This
    should be the person you feel has the most relevant technical background and/or the
    most familiar with the underlying issue or code.
  * The reviewers will be automatically selected based on which module(s) of the code your pull request affects.

The merge request will be reviewed by the core development team before potentially being accepted.

## Setup with Hatch

This project uses [Hatch](https://hatch.pypa.io/latest/).

Although any python environment manager can be used, we recommend using the default eniornment setup by Hatch.

To start using Hatch, it must be installed and accessible from the command line. See the Hatch [installation](https://hatch.pypa.io/latest/install/) for more.

A simple way to install Hatch is to run:

```bash
pip install -g hatch
```

If you can run `hatch -h` then Hatch has been successfully installed.

We recommend setting the `dirs.env` in you hatch config to the following:

```toml
[dirs.env]
virtual = ".hatch"
```

The path to this file can be found by running:

```bash
hatch config find
```

It makes it easier to set the path the the venv in your code editor.

Then run:

```bash
hatch shell
```

This will create the default hatch environment in the project folder. Set your python enviorment path for your editor to `.hatch/eqdsk/bin/python`

## Tests

Run

```bash
hatch run test:tests
```

To run all the tests for the supported python environments
