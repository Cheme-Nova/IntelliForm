# IntelliForm Showcase Pack

This folder captures proof-style IntelliForm runs after the parser/controller reliability upgrade.

- The examples were generated locally on `2026-04-24`.
- They were run through the upgraded controller, not hand-written.
- The test environment did not have a live LLM provider available, so these captures use `LLM_PROVIDER=regex` and the full IntelliForm optimization stack after parsing.

## What Changed Before These Runs

The upgrade applied the most relevant FormulAI learnings to IntelliForm:

- local JSON extraction and repair instead of brittle provider-enforced JSON assumptions
- canonical vertical mapping across parser, API, optimizer, and web UI
- parsed constraints now actually drive optimization
- brief-aware ingredient filtering for obvious contradictions like `silicone-free`, `phosphate-free`, `low-VOC`, and cross-vertical grade leakage
- clearer formulation UX with showcase prompts and parsed-brief visibility

## Included Examples

### 1. OMRI Bioinsecticide Adjuvant

File: [omri-bioinsecticide-adjuvant.json](./omri-bioinsecticide-adjuvant.json)

Why it is worth showing:

- Successful agricultural run with resolved vertical, real blend, cost, eco score, and OMRI-oriented regulatory framing.
- Good example of IntelliForm’s strength in regulatory and sustainability overlays even when the formulation still needs field validation.

What still needs review:

- Tank-mix compatibility and field performance need agronomy validation.
- Some ingredients still require OMRI/manual verification rather than assuming organic acceptance.

### 2. Clean-Label Plant Milk Emulsifier

File: [clean-label-plant-milk-emulsifier.json](./clean-label-plant-milk-emulsifier.json)

Why it is worth showing:

- Successful food-vertical run with explicit allergen flags and a review-required regulatory posture.
- Demonstrates that IntelliForm can return a non-trivial formulation plus useful commercial/regulatory caution instead of a generic “pass.”

What still needs review:

- The current optimizer is better at finding feasible ingredient sets than fully elegant food architecture.
- This should be framed as a first-pass ingredient direction, not a finished commercial beverage system.

### 3. Low-VOC Industrial Degreaser (Infeasible Case)

File: [low-voc-industrial-degreaser-infeasible.json](./low-voc-industrial-degreaser-infeasible.json)

Why it is worth showing:

- This is a strong reliability example: IntelliForm now fails honestly when the filtered industrial pool cannot satisfy the resolved cost, bio-based, and performance constraints.
- It is better proof than hallucinating a fake industrial cleaner.

What still needs review:

- Industrial cleaning still needs deeper formulation architecture and ingredient-package logic in the optimizer.
- This is a good roadmap signal for the next upgrade pass, not a weakness to hide.

## Suggested GitHub Framing

Use language like this in the root README or a launch note:

`These showcase files were captured from the upgraded IntelliForm controller on April 24, 2026. They demonstrate both successful optimization runs and honest infeasibility handling after parser, vertical-mapping, and prompt-reliability upgrades inspired by FormulAI.`
