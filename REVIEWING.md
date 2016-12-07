# Reviewing pull requests

This document provides a high-level, general, and **incomplete** checklist of things that pull request reviewers should be aware of when performing code review. These are guidelines which generally make sense to follow, but they are not intended to be rigid. The checklist mainly consists of things that are specific to the scikit-bio project and that generally apply to incoming pull requests. The checklist is incomplete because it is not possible to describe all things to verify during code review (that depends on what is being reviewed). This document also doesn't attempt to describe *how* to perform code review (there are many online resources for that).

Reviewers are encouraged to keep this document up-to-date as the project evolves and to add anything that's missing.

The checklist is not in any particular order.

## Licensing and attribution

Verify that code being included from external sources has a compatible license and is properly attributed:

- Include the code's license in the top-level `licenses` directory.
- Include a comment with the external code giving attribution and noting the license in the `licenses` directory.
- Any other requirements set by the code's license and/or author.

## Changelog

This is one of the most important points to remember as users will review the changelog to identify changes relevant to them. This is one of the easiest parts to forget in a pull request.

- Note all public (i.e. user-facing) changes in `CHANGELOG.md` under the latest development version section at the top of the file. This includes things like bug fixes, API additions/changes/removal, performance enhancements, etc. The changelog has several subsections for organizing these changes.
- If a corresponding issue exists, it should be linked to from the changelog.
- Use public imports (`skbio.sequence.Sequence` instead of `skbio.sequence._sequence.Sequence`) when documenting import paths in the changelog.
- Internal changes not visible/applicable to users (e.g. refactoring, private methods, etc.) are better suited for commit messages than the changelog.

## Public vs. private API

Be aware of what type of API changes are being made. scikit-bio uses the following conventions to distinguish public vs. private API:

**Public API:** API with an import path that doesn't contain leading underscores in any submodule/subpackage/object names. Users and scikit-bio devs can use public API. Example: `skbio.sequence.Sequence`, `skbio.stats.subsample`

**Package-private API:** API with an import path containing a leading underscore in at least one submodule/subpackage. Users should not import package-private API. Package-private API can be imported anywhere *within* the scikit-bio package. Example: `skbio.util._misc.chunk_str`, `skbio.util._testing.ReallyEqualMixin`

**Private API:** API with object name starting with an underscore. Users should not import private API. Private API should only be imported and used in the module where it is defined. It should not be imported in other parts of the scikit-bio package. Example: `skbio.util._testing._normalize_signs`, `skbio.stats.composition._gram_schmidt_basis`

- Prefer defining private APIs and only promote things to the public API that users need access to.
- Private/package-private API does not need to be decorated with `@experimental`/`@stable`/`@deprecated`, only public API.
- Within scikit-bio, prefer *using* public APIs defined in the package even if private APIs offer the same functionality.

## API stability

Refer to scikit-bio's [API stability docs](http://scikit-bio.org/docs/latest/user/api_stability.html) for the current API lifecycle.

- Prefer making new APIs experimental. This allows for changes to the API without needing to deprecate it first.
- Try to avoid changing stable API if at all possible. If a change is necessary, deprecate the API and remove after 2+ minor releases.
- Under extreme circumstances, a stable API can be changed without deprecation warning to users. This has happened in the past to fix bugs that required an API change. If this happens, discuss with other devs before committing to the change and document this change extensively in the changelog.
- Always document stable/experimental API changes and any deprecated APIs in the changelog.

## Integration/consistency with existing API

When reviewing API changes/additions, look at them in the context of the rest of the codebase and its existing APIs. For example, are there existing parameter names that could be reused for consistency/predictability? Does the new API (or changed API) compose with other relevant APIs? Example: `skbio.stats.distance.anosim` uses a `permutations` parameter, so a new nonparametric function would want to use this name over something like `n` or `num_permutations`.

## Commit messages and merging pull requests

When merging pull requests, use GitHub's "Squash and merge" option to merge a single commit. See [this commit message](https://github.com/biocore/scikit-bio/commit/f3d736aabd717971332781b98d8fde861f354dc3) for an example.

- Rewrite commit message to describe all changes being merged. This usually involves deleting individual commit messages that GitHub includes in the text box.
- Include "fixes #n" text if there's an associated issue to be closed.
- Use [numpy-style commit tags](https://docs.scipy.org/doc/numpy/dev/gitwash/development_workflow.html#writing-the-commit-message) (ENH, BUG, PERF, etc.).
- Include contributors' and reviewers' GitHub usernames in commit message (attribution will be lost on squash).

## Test changes locally

**This step is extremely important.** Pull down the PR changes locally and try out the API as a user would. Try to break it, make sure the docs are complete, etc. Build the docs locally and verify that they render correctly (this is a common pitfall).

## Docs

- Verify the docs follow the instructions in the scikit-bio [documentation guide](https://github.com/biocore/scikit-bio/blob/master/doc/README.md).
- Verify the docs follow [this page](http://scikit-bio.org/docs/latest/development/new_module.html) when adding a new module or subpackage to scikit-bio.
- Public API should have docstrings conforming to [numpydoc standards](https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt). Manual and careful verification of the numpydoc docstrings is currently necessary; they are easy to get wrong and building the docs won't always flag issues. Building the docs and inspecting the rendered output can help with this process.
- Package-private and private APIs do not need to be extensively documented; numpydoc docstrings are not required. Document these APIs as appropriate to help other devs understand the code (code comments are usually better for this anyways).

## CI

- Make sure Travis-CI is passing before merging a pull request.
- Make sure coverage doesn’t drop. Strive to have 100% coverage for new code being merged.

## Unit testing

- Make sure the tests are as complete as possible.
- Check that border cases are tested (e.g. zeros, '', [], None, etc.).
- Check that the base case is tested (`n`), along with the inductive step (`n + 1`).
- Verify that tests cover more than one input data set.
- Make each test case simple, ideally only testing a single thing (follow [Arrange Act Assert](http://wiki.c2.com/?ArrangeActAssert)).
