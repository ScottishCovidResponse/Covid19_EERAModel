# SCRC Software checklist

This checklist is part of ongoing work on a model scoresheet for SCRC models. It relates to software implementation, and assumes that other documents cover questions about model validation, data provenance and quality, quality of science, and policy readiness.

## Software Details

### Model / software name

> Covid19_EERAModel

### Date

> 2020-08-28

### Version identifier

> 0.10.1

## Overall statement

Do we have sufficient confidence in the correctness of the software to trust the results?

This is your overall judgement on the level of confidence based on all the aspects of the checklist. There is no formulaic way to arrive at this overall assessment based on the individual checklist answers but please explain how the measures in place combine to reach this level of confidence and make clear any caveats (eg applies for certain ways of using the software and not others).

> - [ ] Yes
> - [x] Yes, with caveats
> - [ ] No
>
> We are generally quite confident in the correctness of the code. This is based mainly on extensive work on refactoring and reviewing the code by multiple people, helped by regression tests and improvements in the architectural design. But proving correctness via systematic unit and integration testing is still difficult, as that can only really be done via extensive code re-design at a low level to enable effective testing. We are not yet using the data pipeline, but implementation of support for it is ongoing.

## Checklist

Please use a statement from this list: "Sufficiently addressed", "Some work remaining or caveats", or "Needs to be addressed" to begin each response.

Additionally, for each question please explain the situation and include any relevant links (eg tool dashboards, documentation). The sub bullet points are to make the scope of the question clear and should be covered if relevant but do not have to be answered individually.

### Can a run be repeated and reproduce exactly the same results?

- How is stochasticity handled?
- Is sufficient meta-data logged to enable a run to be reproduced: Is the exact code version recorded (and whether the repository was "clean"), including versions of dependent libraries (e.g. an environment.yml file or similar) along with all command line arguments and the content of any configuration files? 
- Is there up-to-date documentation which explains precisely how to run the code to reproduce existing results? 

> - [x] Sufficiently addressed
> - [ ] Some work remaining or caveats
> - [ ] Needs to be addressed
> 
> Project documentation is up to date. The results of a run are a deterministic function of the input data files and the chosen random number seed, and so are reproducible. The software version (Git commit number), and repository clean/dirty flag are logged. Most of the input parameters are logged, but the input files are not copied verbatim to the output directory. Dependency versions in use are logged.

### Are there appropriate tests?  (And are they automated?)

- Are there unit tests? What is covered?
- System and integration tests?  Automated model validation tests?
- Regression tests? (Which show whether changes to the code lead to changes in the output. Changes to the model will be expected to change the output, but many other changes, such as refactoring and adding new features, should not. Having these tests gives confidence that the code hasn't developed bugs due to unintentional changes.)
- Is there CI?
- Is everything you need to run the tests (including documentation) in the repository (or the data pipeline where appropriate)?

> - [ ] Sufficiently addressed
> - [x] Some work remaining or caveats
> - [ ] Needs to be addressed
> 
> We have a set of regression tests which exercise the code in both inference and prediction modes. These have been used mainly as an aid to refactoring, and are not considered confirmation of the code’s correctness. There are a limited set of unit tests which have been developed to test new software modules that we have developed. However, the core model code has effectively no unit tests. As noted above, implementing effective unit tests would require significant re-design of large portions of the core model code. Both our regression tests and unit tests are automated and run in the GitHub CI. Tests can also be run locally on the development host, and there is documentation in the project README on how to do that.

### Are the scientific results of runs robust to different ways of running the code?

- Running on a different machine?
- With different number of processes?
- With different compilers and optimisation levels?
- Running in debug mode?

(We don't require bitwise identical results here, but the broad conclusions after looking at the results of the test case should be the same.) 

> - [x] Sufficiently addressed
> - [ ] Some work remaining or caveats
> - [ ] Needs to be addressed
> 
> The results are robust to running on different environments: we run our regression tests on the CI both in Linux and MacOS environments. We have not run the code using different compiler configurations (e.g. debug versus optimised), but we do not believe that this would menaingfully change the results. Our code is single-threaded at the moment, and doesn’t use any machine-specific features, so there should be little potential for variation across different machines.

### Has any sort of automated code checking been applied?

- For C++, this might just be the compiler output when run with "all warnings". It could also be more extensive static analysis. For other languages, it could be e.g. pylint, StaticLint.jl, etc.
- If there are possible issues reported by such a tool, have they all been either fixed or understood to not be important?

> - [x] Sufficiently addressed
> - [ ] Some work remaining or caveats
> - [ ] Needs to be addressed
> 
> We switch on compiler warnings and use a variety of static analysis tools (cppcheck, clang-tidy and a few others). These are all integrated with the GitHUb CI, amd all of the output is visible on the project dashboard (https://scottishcovidresponse.github.io/Covid19_EERAModel/). A pass through the code has been done to eliminate most of the warnings: the remaining warnings are largely concerning cosmetic issues and are considered a low priority for resolving. 

### Is the code clean, generally understandable and readable and written according to good software engineering principles?

- Is it modular?  Are the internal implementation details of one module hidden from other modules?
- Commented where necessary?
- Avoiding red flags such as very long functions, global variables, copy and pasted code, etc.?

> - [x] Sufficiently addressed
> - [ ] Some work remaining or caveats
> - [ ] Needs to be addressed
> 
> The project sources are in a good state overall. After extensive refactoring, there is good use of modularity, inline documentation, well-defined interfaces hiding implementation details, and fairly extensive use of modern C++ idioms. There is some code duplication specifically in the implementations of the three epidemiological models: this is deliberate, to allow maximum flexibility in modifying individual models independently of each other. 

### Is there sufficient documentation?

- Is there a readme?
- Does the code have user documentation?
- Does the code have developer documentation?
- Does the code have algorithm documentation? e.g. something that describes how the model is actually simulated, or inference is performed?
- Is all the documentation up to date? 

> - [x] Sufficiently addressed
> - [ ] Some work remaining or caveats
> - [ ] Needs to be addressed
> 
> Documentation exists at all of these levels. There is high-level documentation (in the `docs/` directory) describing how the model works from an epidemiological point of view. There is user documentation (`README.md`) explaining how to build and run the code, and what inputs are required, and what outputs are produced. There is also inline code documentation which is exported via Doxygen, and can be viewed on the project dashboard (https://scottishcovidresponse.github.io/Covid19_EERAModel/site/doxygen-docs.html).

### Is there suitable collaboration infrastructure?

- Is the code in a version-controlled repository?
- Is there a license?
- Is an issue tracker used?
- Are there contribution guidelines?

> - [x] Sufficiently addressed
> - [ ] Some work remaining or caveats
> - [ ] Needs to be addressed
> 
> Project code is stored on GitHub at https://github.com/ScottishCovidResponse/Covid19_EERAModel. The code
is licensed under the BSD-2-Clause license. Issues are recorded on the SCRC issue tracker (https://github.com/ScottishCovidResponse/SCRCIssueTracking). Contribution guidelines are provided in `CONTRIBUTING.md`.

### Are software dependencies listed and of appropriate quality?

> - [x] Sufficiently addressed
> - [ ] Some work remaining or caveats
> - [ ] Needs to be addressed
> 
> Libraries are listed in the project README and are of appropriate quality. The principal dependency is the GNU Scientific Library (http://www.gnu.org/software/gsl/), which is a high-quality numerical library for C/C++. Other dependencies are listed in the project README.

### Is input and output data handled carefully?

- Does the code use the data pipeline for all inputs and outputs?
- Is the code appropriately parameterized (i.e. have hard coded parameters been removed)?

> - [ ] Sufficiently addressed
> - [x] Some work remaining or caveats
> - [ ] Needs to be addressed
> 
> The code is not yet using the data pipeline. We have run a successful build against the pipeline API, but this is purely a proof of principle demonstration. Trial EERA model parameters have been uploaded to the data pipeline repositories. The current focus is on re-design of the I/O portions of the model code to support making API calls to read and write the appropriate data items. The code is entirely parameterised: there are no hard-coded parameters.
