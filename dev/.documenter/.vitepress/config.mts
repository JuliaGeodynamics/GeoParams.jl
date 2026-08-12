import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import path from 'path'

// console.log(process.env)

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/GeoParams.jl/dev/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [{ text: 'Home', link: '/index' },
{ text: 'User Guide',  items: [
{text: 'GeoUnit', link: '/man/geounit'},
{text: 'Nondimensionalization', link: '/man/nondimensionalize'},
{text: 'Material Parameters',  items: [
{text: 'Overview', link: '/man/materialparameters'},
{text: 'Permeability', link: '/man/permeability'},
{text: 'Heat Capacity', link: '/man/heatcapacity'},
{text: 'Conductivity', link: '/man/conductivity'},
{text: 'Latent heat', link: '/man/latentheat'},
{text: 'Radioactive heat', link: '/man/radioactiveheating'},
{text: 'Shear heating', link: '/man/shearheating'},
{text: 'Gravity', link: '/man/gravity'},
{text: 'Partial Melting', link: '/man/melting'},
{text: 'Density', link: '/man/density'},
{text: 'Solubility', link: '/man/solubility'}
]},
{text: 'Constitutive Relationships',  items: [
{text: 'Creep laws', link: '/man/creeplaws'},
{text: 'Custom rheology', link: '/man/customrheology'},
{text: 'Viscosity', link: '/man/viscosity'},
{text: 'Elasticity', link: '/man/elasticity'},
{text: 'Plasticity', link: '/man/plasticity'}
]},
{text: 'Chemical Diffusion',  items: [
{text: 'Computational routines', link: '/man/chemicaldiffusion'},
{text: 'Garnet', link: '/man/Garnet'},
{text: 'Melt', link: '/man/Melt'},
{text: 'Olivine', link: '/man/Olivine'},
{text: 'Rutile', link: '/man/Rutile'}
]},
{text: 'TAS classification', link: '/man/TASclassification'},
{text: 'Zircon Ages', link: '/man/zirconages'},
{text: 'Phase Diagrams', link: '/man/phasediagrams'},
{text: 'Seismic Velocity', link: '/man/seismicvelocity'},
{text: '1D Strength Envelope', link: '/man/strengthenvelope'}
] },
{ text: 'Plotting', link: '/man/plotting' },
{ text: 'List of functions', link: '/man/listfunctions' },
{ text: 'Contributing', link: '/man/contributing' }]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker',
  }
]
// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/GeoParams.jl/dev/', // TODO: replace this in makedocs!
  title: 'GeoParams.jl',
  description: 'Documentation for GeoParams.jl',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../1', // This is required for MarkdownVitepress to work correctly...

  head: [
    ['meta', { name: 'robots', content: 'noindex, nofollow' }],
    ['link', { rel: 'icon', href: `${baseTemp.base}favicon.ico` }],
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  vite: {
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    build: {
      assetsInlineLimit: 0, // so we can tell whether we have created inlined images or not, we don't let vite inline them
    },
    optimizeDeps: {
      exclude: [
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ],
    },
    ssr: {
      noExternal: [
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ],
    },
    assetsInclude: ['**/*.PNG'], // Add this line to handle .PNG files as assets
  },

  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"
    },
  },
  themeConfig: {
    outline: 'deep',
    // https://vitepress.dev/reference/default-theme-config
    logo: { src: '/logo.png', width: 24, height: 24},
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [{ text: 'Home', link: '/index' },
{ text: 'User Guide', collapsed: false, items: [
{text: 'GeoUnit', link: '/man/geounit'},
{text: 'Nondimensionalization', link: '/man/nondimensionalize'},
{text: 'Material Parameters', collapsed: false, items: [
{text: 'Overview', link: '/man/materialparameters'},
{text: 'Permeability', link: '/man/permeability'},
{text: 'Heat Capacity', link: '/man/heatcapacity'},
{text: 'Conductivity', link: '/man/conductivity'},
{text: 'Latent heat', link: '/man/latentheat'},
{text: 'Radioactive heat', link: '/man/radioactiveheating'},
{text: 'Shear heating', link: '/man/shearheating'},
{text: 'Gravity', link: '/man/gravity'},
{text: 'Partial Melting', link: '/man/melting'},
{text: 'Density', link: '/man/density'},
{text: 'Solubility', link: '/man/solubility'}
]},
{text: 'Constitutive Relationships', collapsed: false, items: [
{text: 'Creep laws', link: '/man/creeplaws'},
{text: 'Custom rheology', link: '/man/customrheology'},
{text: 'Viscosity', link: '/man/viscosity'},
{text: 'Elasticity', link: '/man/elasticity'},
{text: 'Plasticity', link: '/man/plasticity'}
]},
{text: 'Chemical Diffusion', collapsed: false, items: [
{text: 'Computational routines', link: '/man/chemicaldiffusion'},
{text: 'Garnet', link: '/man/Garnet'},
{text: 'Melt', link: '/man/Melt'},
{text: 'Olivine', link: '/man/Olivine'},
{text: 'Rutile', link: '/man/Rutile'}
]},
{text: 'TAS classification', link: '/man/TASclassification'},
{text: 'Zircon Ages', link: '/man/zirconages'},
{text: 'Phase Diagrams', link: '/man/phasediagrams'},
{text: 'Seismic Velocity', link: '/man/seismicvelocity'},
{text: '1D Strength Envelope', link: '/man/strengthenvelope'}
] },
{ text: 'Plotting', link: '/man/plotting' },
{ text: 'List of functions', link: '/man/listfunctions' },
{ text: 'Contributing', link: '/man/contributing' }]
,
    editLink: { pattern: "https://github.com/JuliaGeodynamics/GeoParams.jl/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/JuliaGeodynamics/GeoParams.jl' },
      { icon: 'slack', link: 'https://julialang.org/slack/' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()} ⋅ The GeoParams Development Team.`
    }
  }
})
