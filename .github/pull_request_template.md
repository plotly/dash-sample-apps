Issue for app: #[issue number here]

# App pull request
- [ ] This is a new app
- [ ] I am improving an existing app (redesigns/code "makeovers")

## About

- Playground deployment URL (new version):
- Current gallery app URL: (delete this line if inapplicable)
- Python app repository link: (delete this line if you are working on a Python app)

## Workflow

- [ ] I have created a branch in the appropriate monorepo, and the
  elements necessary for successful deployment are in place.
- [ ] If the app is a redesigned and/or restyled version of an
  existing gallery app, I've summarized the changes requested in the
  appropriate Streambed issue and confirm that they have been applied.
- [ ] If the app is on the Dash Gallery portal, I have added a link to
  the GitHub repository for the source code in the portal description.
- [ ] If the app is a reimplementation of a Python gallery app for the
  DashR gallery, the app in this PR mimics, as closely as possible,
  the style and functionality of the existing app.=
- [ ] I have removed *all Google Analytics code* from the app's
  `assets/` folder.

## The pre-review review

I have addressed all of the following questions:

- [ ] Does everything in my code serve some purpose? (I have removed
  any dead and/or irrelevant code.)
- [ ] Does everything in my code have a clear purpose? (My code is
  readable and, where it isn't, it has been commented appropriately.)]
- [ ] Am I reinventing the wheel? (I have used appropriate packages to
  lessen the volume of code that needs to be maintained.)

## Post PR (at merge time)

- [ ] When you are merging, make sure to write one of the following tags in the commit message (or it will default to patch):
  - `#patch` - An app has been updated or fixed
  - `#minor` - A new app has been added, or an app has been significantly reworked
  - `#major` - Breaking changes, make sure to discuss with dash-core before using this tag