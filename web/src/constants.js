export const VERTICAL_OPTIONS = [
  { value: 'personal_care', label: 'Personal Care' },
  { value: 'industrial', label: 'Industrial' },
  { value: 'agricultural', label: 'Agricultural' },
  { value: 'pharmaceutical', label: 'Pharmaceutical' },
  { value: 'food', label: 'Food & Beverage' },
  { value: 'fabric_laundry', label: 'Fabric & Laundry' },
  { value: 'paint_coatings', label: 'Paint & Coatings' },
]

export const PUBLIC_VERTICAL_GUIDES = {
  personal_care: {
    status: 'beta',
    label: 'Beta',
    message: 'Personal care runs in the public edition and usually returns a feasible first-pass blend, but the architecture still needs refinement before we would treat it as showcase-ready for premium cosmetic claims.',
    prompts: [
      {
        title: 'Clear Mild Body Wash',
        vertical: 'personal_care',
        text: 'Sulfate-free clear body wash with mild foam and rinse-clean feel under $5/kg.',
        note: 'Currently the cleanest public personal-care starter. Use as a beta demo rather than a final cosmetic claim set.',
      },
    ],
  },
  industrial: {
    status: 'beta',
    label: 'Beta',
    message: 'Industrial remains a beta public vertical. This starter is included because it fails honestly under the current public ingredient set, which is still better than bluffing a fake industrial blend.',
    prompts: [
      {
        title: 'Low-Odor Hard-Surface Cleaner',
        vertical: 'industrial',
        text: 'General-purpose hard-surface cleaner concentrate for industrial maintenance, low odor, under $6/kg.',
        note: 'Currently a parser-and-feasibility test prompt. Expect an honest infeasible result more often than a polished final blend.',
      },
    ],
  },
  agricultural: {
    status: 'validated',
    label: 'Validated',
    message: 'These prompts are the strongest current public starters for showing IntelliForm’s parser, optimization, and reporting flow in a way that usually lands cleanly.',
    prompts: [
      {
        title: 'OMRI Bioinsecticide Adjuvant',
        vertical: 'agricultural',
        text: 'OMRI-listed bioinsecticide adjuvant for organic crops, bee-safe, under $8/kg, improves spray coverage.',
        note: 'Strong current agricultural demo with a clear adjuvant use case.',
      },
      {
        title: 'Bio-Based Biostimulant Blend',
        vertical: 'agricultural',
        text: 'Biostimulant blend for enhanced root growth, 100% bio-based, soil-safe, and minimal persistence.',
        note: 'Good for showing sustainability posture and a successful optimization outcome.',
      },
    ],
  },
  pharmaceutical: {
    status: 'beta',
    label: 'Beta',
    message: 'Pharmaceutical remains a beta public vertical. The starter below is kept for transparent testing, but the public excipient architecture is not yet ready for claim-heavy showcase use.',
    prompts: [
      {
        title: 'Tablet Film Coating System',
        vertical: 'pharmaceutical',
        text: 'Immediate-release tablet film coating system, water-based, pilot-friendly.',
        note: 'Useful for testing parser and feasibility behavior, but not yet a strong public success case.',
      },
    ],
  },
  food: {
    status: 'validated',
    label: 'Validated',
    message: 'Food outputs are not yet claim-ready, but these prompts consistently demonstrate the current public engine’s ability to parse and return feasible first-pass concepts.',
    prompts: [
      {
        title: 'Plant-Milk Emulsifier Blend',
        vertical: 'food',
        text: 'Clean-label emulsifier blend for plant-based milk, GRAS approved, under $4/kg, neutral taste, stable at pasteurization temperature.',
        note: 'Useful for a clean-label, functional-ingredient style demo.',
      },
      {
        title: 'Natural Bakery Preservation System',
        vertical: 'food',
        text: 'Natural preservative system for bakery applications with clean-label positioning and mold control under $4/kg.',
        note: 'Shows constraint handling in a recognizable food formulation brief.',
      },
    ],
  },
  fabric_laundry: {
    status: 'validated',
    label: 'Validated',
    message: 'Fabric and laundry is currently one of the better-performing public verticals for a first interactive demo.',
    prompts: [
      {
        title: 'Cold-Water Laundry Detergent',
        vertical: 'fabric_laundry',
        text: 'Phosphate-free laundry detergent for cold-water washing with stain focus and improved sustainability profile under $4/kg.',
        note: 'Good current balance of feasibility and category relevance.',
      },
      {
        title: 'Bio-Focused Laundry Liquid',
        vertical: 'fabric_laundry',
        text: 'Concentrated liquid detergent, phosphate-free, cold-water active, stain-focused, and at least 85% bio-based.',
        note: 'Stronger sustainability target for users who want to stress-test the public optimizer.',
      },
    ],
  },
  paint_coatings: {
    status: 'beta',
    label: 'Beta',
    message: 'Paint and coatings is still a beta public vertical. This starter is included for transparent testing while the architecture and ingredient coverage are tuned further.',
    prompts: [
      {
        title: 'Waterborne Coating Aid',
        vertical: 'paint_coatings',
        text: 'Waterborne coating aid for interior paint with low odor and easy application.',
        note: 'Currently best used as a public beta feasibility check rather than a proof-quality optimization case.',
      },
    ],
  },
}
