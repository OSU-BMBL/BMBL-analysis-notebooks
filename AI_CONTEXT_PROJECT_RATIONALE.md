# RATIONALE: Workflow-Specific AI Context Enhancement

**Project**: BMBL AI-Friendly Documentation - Phase 5  
**Status**: Implementation Phase  
**Started**: April 2026  
**Last Updated**: April 2026

---

## Executive Summary

This project creates deep contextual understanding for AI assistants working with specific BMBL workflows. Following the successful completion of Phases 1-4 (documentation, dependencies, CI, reproducibility), Phase 5 focuses on enabling AI to provide accurate, workflow-specific assistance rather than generic suggestions.

---

## Problem Statement

### Current Gap
While Phases 1-4 provided:
- ✅ High-level repo navigation (AGENTS.md)
- ✅ Design philosophy (RATIONALE.md)
- ✅ Automated validation (CI)
- ✅ Reproducible environments (Docker)

**What's missing**: AI lacks **deep, workflow-specific context** to provide accurate help when users work on specific analyses.

### Real-World Impact
**Scenario**: User asks AI to "help modify the scRNAseq workflow to change normalization"

**Without context**:
- AI suggests changing wrong variable
- Misses downstream dependencies
- Doesn't know current method or why it was chosen
- May suggest incompatible approach

**With context** (what we're building):
- AI knows exact file and line
- Understands impact on subsequent steps
- Suggests compatible alternatives
- Provides testing guidance

---

## Solution: `.ai_context.md` Files

### Concept
Create a context file in each major workflow directory that provides AI assistants with:
1. **Data flow understanding** - How data transforms step-by-step
2. **Common modification patterns** - Typical changes and implementations
3. **Gotchas and warnings** - Things that commonly break
4. **File relationships** - Dependencies between scripts
5. **Testing guidance** - How to verify changes work

### File Location
```
workflow_name/
├── README.md              # Human documentation
├── .ai_context.md         # AI deep context ← NEW
├── 0_install_packages.R
├── 1_preprocess.rmd
└── ...
```

---

## Why This Approach Won

### Comparison Matrix

| Approach | Deep Context | Maintainable | Scalable | Human-Readable | Our Choice |
|----------|--------------|--------------|----------|----------------|------------|
| `.ai_context.md` | ✅ Yes | ✅ Yes | ✅ Yes | ✅ Yes | ✅ **WINNER** |
| Inline comments | ✅ Yes | ❌ Clutters code | ❌ Hard | ✅ Yes | ❌ |
| Prompt templates | ❌ No | ✅ Yes | ✅ Yes | ✅ Yes | ❌ |
| .cursorrules | ❌ Generic | ✅ Yes | ✅ Yes | ✅ Yes | ❌ |
| JSON metadata | ✅ Yes | ⚠️ OK | ✅ Yes | ❌ No | ❌ |

### Key Advantages
1. **Just-in-time context**: AI reads only when working on that workflow
2. **Markdown format**: Easy to write, maintain, version control
3. **Human + AI readable**: Lab members can review and update
4. **Proven pattern**: Builds on success of AGENTS.md
5. **Scalable**: Start with 5 workflows, expand gradually

---

## Implementation Plan

### Phase 1: Template Creation & Testing (Current)
**Goal**: Create template and validate on one workflow

**Tasks**:
1. [ ] Create `.ai_context_TEMPLATE.md` in repository root
2. [ ] Apply template to `scRNAseq_general_workflow/`
3. [ ] Test with actual AI queries
4. [ ] Refine template based on learnings

**Timeline**: 1 session

**Success Criteria**:
- Template covers all necessary sections
- AI can answer workflow-specific questions accurately
- Format is maintainable

---

### Phase 2: Common Recipes (After Phase 1 validation)
**Goal**: Create reusable recipes for common tasks

**Structure**:
```
_common/
├── functions.R              # Existing
└── ai_recipes.md           # NEW - General recipes
    
[workflow]/
├── .ai_context.md          # References recipes + adds specifics
```

**Recipe Categories**:
- Data manipulation (filtering, subsetting)
- Visualization (adding plots, changing aesthetics)
- Export (to different formats)
- Debugging (memory, errors)
- Parameters (changing analysis settings)

**Timeline**: 1 session (after Phase 1)

---

### Phase 3: Roll Out to Remaining 4 Workflows
**Goal**: Apply refined template to all target workflows

**Target Workflows**:
1. ✅ `scRNAseq_general_workflow/` (done in Phase 1)
2. [ ] `scRNAseq_trajectory_Slingshot/`
3. [ ] `scATACseq_general_workflow/`
4. [ ] `RNAseq_nfcore_workflow/`
5. [ ] `ST_general_workflow/`

**Timeline**: 1-2 sessions

---

### Phase 4: Integration & Documentation
**Goal**: Make AI context files discoverable and maintainable

**Tasks**:
- [ ] Update AGENTS.md to mention .ai_context.md
- [ ] Add validation check in validate_repo.R
- [ ] Create contributor guide for adding context files
- [ ] Update RATIONALE.md with Phase 5 documentation

**Timeline**: 1 session

---

## Template Structure

### Required Sections

```markdown
# AI Context: [Workflow Name]

## Quick Summary
- One paragraph overview
- Primary analysis type
- Typical runtime and resources

## Data Flow
Step-by-step with file locations:
1. Raw input → 2. QC → 3. Normalization → 4. Analysis → 5. Output

## Common Modifications
| Task | File | Line/Variable | Notes |
|------|------|---------------|-------|
| Change resolution | 2_cluster.rmd | `resolution = 0.8` | Affects cluster count |

## Gotchas & Warnings
- Must run in order (1 → 2 → 3)
- Common error X means Y
- Don't change Z without updating W

## File Relationships
- File A outputs → File B inputs
- Shared variables between files

## Testing
- Quick test: use `pbmc_small` dataset
- Full test: takes X min, Y GB RAM
- Expected outputs: ...

## External Dependencies
- OSC modules required
- Reference databases
- External tools

## See Also
- Link to general recipes in _common/
- Related workflows
```

---

## Target Workflows Detail

### 1. scRNAseq_general_workflow (Phase 1)
**Why first**: Most popular, well-understood, good test case
**Key aspects**: Seurat, clustering, annotation
**Expected recipes**: Normalization, clustering, markers, visualization

### 2. scRNAseq_trajectory_Slingshot
**Why**: Different analysis type (pseudotime vs static)
**Key aspects**: Trajectory inference, gene dynamics
**Expected recipes**: Lineage selection, gene trends, branch analysis

### 3. scATACseq_general_workflow
**Why**: Different data modality (ATAC vs RNA)
**Key aspects**: Peak calling, motif analysis, integration
**Expected recipes**: Peak thresholds, motif databases, multiome

### 4. RNAseq_nfcore_workflow
**Why**: Different approach (bulk vs single-cell)
**Key aspects**: DESeq2, pathways, nfcore pipeline
**Expected recipes**: Contrast design, pathway analysis, batch effects

### 5. ST_general_workflow
**Why**: Emerging technology (spatial)
**Key aspects**: Spatial coordinates, neighborhood analysis
**Expected recipes**: Spot selection, spatial features, image integration

---

## Success Metrics

### Quantitative
- [ ] AI can answer 80%+ of workflow-specific questions correctly
- [ ] Time to resolve workflow questions: < 5 minutes
- [ ] Reduction in AI suggestion errors: 50%+

### Qualitative
- [ ] Users report AI is "more helpful" for workflow tasks
- [ ] New lab members can understand workflows faster
- [ ] Fewer back-and-forth clarifications needed

---

## Maintenance Strategy

### When to Update
- Workflow structure changes
- New common gotcha discovered
- External dependency added/removed
- Quarterly review (calendar reminder)

### Who Maintains
- **Primary**: Original workflow author
- **Secondary**: Anyone modifying the workflow
- **Community**: Based on AI interaction feedback

### How to Update
1. Edit `.ai_context.md` in workflow directory
2. Test changes with actual AI query
3. Commit with message: "Update AI context for [workflow]"
4. Update if validation check fails

---

## Integration with Existing Infrastructure

### CI/CD (Phase 3)
Add check in `validate_repo.R`:
```r
# Warn if workflow has R code but no .ai_context.md
if (has_r_code && !has_ai_context) {
  warning("Consider adding .ai_context.md for better AI assistance")
}
```

### Documentation (Phase 1)
Update `AGENTS.md`:
```markdown
## Working with Workflows
Each workflow has an `.ai_context.md` file with AI-specific guidance
on data flow, common modifications, and troubleshooting.
```

### Validation (Phase 3)
Check `.ai_context.md` is valid markdown and has required sections.

---

## Future Expansion

### Beyond 5 Workflows
- Roll out to remaining 42+ workflows gradually
- Priority: most-used workflows first
- Can be done incrementally over time

### Advanced Features (Future Phases)
- **Structured metadata**: JSON schema for machine parsing
- **Interactive queries**: AI can ask clarifying questions
- **Version tracking**: Track how workflow changes affect context
- **Auto-generation**: Extract some info from code automatically

---

## Questions & Notes

### Open Questions
- Should we include code snippets in context files?
- How detailed should "common modifications" be?
- Should we track which recipes are most-used?

### Known Limitations
- AI must be instructed to read .ai_context.md
- Won't help if user asks about workflow without specifying which one
- Requires maintenance when workflows change

---

## References

- AGENTS.md - High-level repo navigation
- RATIONALE.md - Design philosophy (all phases)
- validate_repo.R - Validation infrastructure
- Workflow READMEs - Human documentation

---

**Next Step**: Create template and apply to scRNAseq_general_workflow

**Status**: Ready to implement Phase 1
