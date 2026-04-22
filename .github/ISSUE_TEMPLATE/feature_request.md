---
name: Feature Request
description: Suggest a new feature or enhancement
title: "[FEATURE] "
labels: ["enhancement", "triage"]
body:
  - type: markdown
    attributes:
      value: |
        Thanks for suggesting a feature! Please describe your idea below.

  - type: textarea
    id: problem
    attributes:
      label: Problem/Need
      description: What problem does this solve? What need does it address?
      placeholder: I'm always frustrated when... / I need to be able to...
    validations:
      required: true

  - type: textarea
    id: solution
    attributes:
      label: Proposed Solution
      description: Describe the feature you'd like to see
      placeholder: A clear description of what you want to happen...
    validations:
      required: true

  - type: textarea
    id: alternatives
    attributes:
      label: Alternatives Considered
      description: What alternative solutions have you considered?
      placeholder: Other approaches you've tried or thought about...

  - type: dropdown
    id: scope
    attributes:
      label: Scope
      description: Where would this feature fit?
      options:
        - New workflow (entirely new analysis type)
        - Enhancement to existing workflow
        - Utility/helper function
        - Documentation improvement
        - Other
    validations:
      required: true

  - type: textarea
    id: context
    attributes:
      label: Additional Context
      description: Any other information, mockups, or examples
