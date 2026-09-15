# AI-Assisted Prime Factorisation and Number Theory Research

This repository contains an ongoing AI-assisted research project exploring **prime factorisation, number theory, integer identities, and related mathematical structures**.

The goal is to use AI as a research assistant for exploring mathematical ideas, generating hypotheses, designing computational experiments, and helping formalize results.

The project is **experimental**: computational evidence is used to discover and investigate patterns, but experimental results are not considered proofs by themselves.

## How it works

The research follows an iterative process:

```text
Research papers / mathematical knowledge
                ↓
        Understand the area
                ↓
        Generate hypotheses
                ↓
       Prioritize experiments
                ↓
       Implement experiments
                ↓
          Gather data
                ↓
     Search for identities
                ↓
      Form conjectures
                ↓
       Prove or refute them
                ↓
       Document the results
                ↓
          Write papers
```

### Research workflow

1. Feed the AI relevant mathematical knowledge, such as research papers, books, notes, and previously discovered results.

2. Ask the AI to explain the relevant mathematical area, definitions, constructions, and identities.

3. Ask the AI to propose possible mathematical directions and experiments.

4. Create an experiment priority list based on mathematical interest, expected usefulness, and computational feasibility.

5. Ask the AI to help implement a Python experiment for testing a specific hypothesis.

6. Run the experiment and collect the generated results.

7. Analyze the results and look for mathematical patterns, identities, counterexamples, or possible conjectures.

8. Investigate promising identities mathematically and attempt to prove or disprove them.

9. Repeat the experimental process until the mathematical structure is sufficiently understood.

10. After a substantial group of experiments, consolidate the results into a LaTeX research paper.

11. When a research thread becomes sufficiently large, start a new AI conversation and continue from the documented results rather than relying on the previous conversation history.

## Experiments

The `experiments/` directory contains the computational experiments used during the research.

Each experiment should ideally document:

* the experiment number;
* the mathematical hypothesis being tested;
* the construction or formula being investigated;
* how the test data is generated;
* the number and size of test cases;
* the experimental results;
* unexpected observations;
* failures and counterexamples;
* conclusions;
* and, where applicable, the mathematical identity or conjecture discovered.

Experiments are intended to be **reproducible** whenever possible.

A computational result is treated as evidence for a hypothesis, not as a mathematical proof.

## Conversations

The `conversations/` directory contains the AI conversations that led to experiments, mathematical ideas, conjectures, and papers.

These conversations provide the research history and show how ideas developed over time.

They can contain:

* mathematical discussions;
* explanations of existing research;
* experiment proposals;
* experiment analysis;
* failed approaches;
* discovered identities;
* proof attempts;
* and research directions.

Failed experiments and abandoned approaches are intentionally preserved when they are useful for understanding the research process.

## Papers

The `papers/` directory contains the mathematical papers produced during the project.

Papers may contain results originating from multiple experiments.

A typical progression is:

```text
Experiment
    ↓
Observation
    ↓
Conjecture
    ↓
Mathematical investigation
    ↓
Proof / counterexample
    ↓
Research result
    ↓
Paper
```

Only results that survive mathematical scrutiny should be presented as established mathematical results.

## Philosophy

The purpose of using AI in this project is not to assume that the AI is correct.

Instead, AI is used as a tool for:

* exploring large mathematical search spaces;
* finding potentially interesting connections;
* generating computational experiments;
* suggesting alternative formulations;
* identifying possible patterns;
* assisting with symbolic manipulation;
* and helping organize mathematical research.

AI-generated claims must therefore be independently checked.

A useful distinction throughout the repository is:

**Observation → Conjecture → Proof**

An observation obtained from computation is not automatically a conjecture, and a conjecture is not automatically a theorem.

## Reproducibility

Experiments should generate their own test data rather than depending on previous experiment output whenever practical.

Experiment scripts should clearly identify themselves when executed, for example:

```text
START EXPERIMENT 23

...

FINISHED EXPERIMENT 23
```

This makes long-running experiment logs easier to identify and compare.

## Research history

The repository intentionally preserves the development of the research rather than only its final results.

This includes:

* successful experiments;
* failed experiments;
* incorrect conjectures;
* counterexamples;
* discarded approaches;
* intermediate identities;
* and earlier versions of mathematical ideas.

The failed paths are part of the research record and can sometimes provide useful information for future experiments.

# Credits

* **OpenAI — ChatGPT**: AI research assistance.
* The mathematical literature and researchers whose published work is used as input for the research.

## Disclaimer

This repository contains exploratory mathematical research.

AI-generated explanations, conjectures, code, and mathematical statements may contain errors. Computational experiments should be independently reproduced, and mathematical claims should be verified through rigorous reasoning or proof before being treated as established results.
