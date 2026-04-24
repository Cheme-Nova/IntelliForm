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
    message: 'Personal care is available in the public edition, but the current engine still needs stronger conditioner and cleanser architecture logic before we recommend it as a first-click demo vertical.',
    prompts: [],
  },
  industrial: {
    status: 'beta',
    label: 'Beta',
    message: 'Industrial prompts remain available, but the public ingredient set still struggles with realistic low-VOC and release-agent targets. We do not recommend industrial as a first demo click yet.',
    prompts: [],
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
    message: 'Pharmaceutical logic is present, but the public edition still needs more trustworthy excipient architecture before it should be showcased as a first-use vertical.',
    prompts: [],
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
    message: 'Paint and coatings remains a beta vertical in the public edition. We recommend testing it later, after the public formulation architecture is tuned further.',
    prompts: [],
  },
}
