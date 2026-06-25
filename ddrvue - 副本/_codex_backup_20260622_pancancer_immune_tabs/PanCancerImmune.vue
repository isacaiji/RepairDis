<template>
  <div class="pancancer-immune">
    <section class="hero-panel">
      <div>
        <p class="eyebrow">MO-DDRscore-derived pan-cancer immune landscape</p>
        <h2>Pan-cancer DDR-state immune analysis</h2>
        <p class="hero-text">
          RepairDis provides precomputed immune analyses based on MO-DDRscore high- and low-state
          stratification across TCGA cancer types. Results cover tumor microenvironment scores,
          immune cell remodeling, checkpoint biomarkers and TIDE-predicted immunotherapy response.
        </p>
      </div>
      <div class="summary-grid">
        <div class="summary-card">
          <strong>33</strong>
          <span>cancer types</span>
        </div>
        <div class="summary-card">
          <strong>4</strong>
          <span>immune views</span>
        </div>
        <div class="summary-card">
          <strong>Static</strong>
          <span>precomputed results</span>
        </div>
      </div>
    </section>

    <section class="figure-grid">
      <article v-for="item in figureItems" :key="item.title" class="figure-card">
        <div class="figure-card-header">
          <div>
            <p class="panel-tag">{{ item.tag }}</p>
            <h3>{{ item.title }}</h3>
          </div>
          <div class="action-row">
            <a :href="item.pdf" target="_blank" rel="noopener">Open PDF</a>
            <a :href="item.table" download>Data</a>
          </div>
        </div>
        <p class="figure-desc">{{ item.description }}</p>
        <div class="figure-preview">
          <img :src="item.image" :alt="item.title" loading="lazy" />
        </div>
      </article>
    </section>

    <section class="download-panel">
      <div>
        <p class="eyebrow">Download</p>
        <h3>Reusable precomputed immune result tables</h3>
        <p>
          These files are provided as analysis-ready tables and can be reused for external plotting,
          validation, or manuscript supplement preparation.
        </p>
      </div>
      <div class="download-list">
        <a v-for="file in tableItems" :key="file.href" :href="file.href" download>
          {{ file.label }}
        </a>
      </div>
    </section>
  </div>
</template>

<script setup>
const base = '/repairdis/pancancer-immune'

const figureItems = [
  {
    tag: 'TME score',
    title: 'Global tumor microenvironment',
    description: 'ESTIMATE-derived tumor purity, stromal score, immune score and ESTIMATE score across MO-DDRscore states.',
    image: `${base}/figures/Figure2A_ESTIMATE_landscape_locked.png`,
    pdf: `${base}/figures/Figure2A_ESTIMATE_landscape_locked.pdf`,
    table: `${base}/tables/Figure2A_ESTIMATE_locked_plot_data.csv`
  },
  {
    tag: 'Cell type',
    title: 'Immune cell-type remodeling',
    description: 'Representative immune cell scores summarized between MO-DDRscore-low and -high tumors across cancer types.',
    image: `${base}/figures/Figure2B_TME_cell_dumbbell_locked.png`,
    pdf: `${base}/figures/Figure2B_TME_cell_dumbbell_locked.pdf`,
    table: `${base}/tables/Figure2B_immune_cell_locked_plot_data.csv`
  },
  {
    tag: 'Checkpoint',
    title: 'Checkpoint biomarker landscape',
    description: 'Checkpoint gene expression patterns associated with MO-DDRscore high-low contrast across pan-cancer cohorts.',
    image: `${base}/figures/Figure2D_checkpoint_landscape_locked.png`,
    pdf: `${base}/figures/Figure2D_checkpoint_landscape_locked.pdf`,
    table: `${base}/tables/Figure2D_checkpoint_locked_plot_data.csv`
  },
  {
    tag: 'TIDE',
    title: 'Predicted immunotherapy response',
    description: 'TIDE-predicted responder proportions in MO-DDRscore-low and -high tumors across cancer types.',
    image: `${base}/figures/Figure2E_TIDE_response_red_bottom_locked.png`,
    pdf: `${base}/figures/Figure2E_TIDE_response_red_bottom_locked.pdf`,
    table: `${base}/tables/Figure2E_TIDE_response_red_bottom_plot_data.csv`
  }
]

const tableItems = [
  {
    label: 'ESTIMATE landscape table',
    href: `${base}/tables/Figure2A_ESTIMATE_locked_plot_data.csv`
  },
  {
    label: 'Immune cell remodeling table',
    href: `${base}/tables/Figure2B_immune_cell_locked_plot_data.csv`
  },
  {
    label: 'Checkpoint biomarker table',
    href: `${base}/tables/Figure2D_checkpoint_locked_plot_data.csv`
  },
  {
    label: 'TIDE response proportion table',
    href: `${base}/tables/Figure2E_TIDE_response_red_bottom_plot_data.csv`
  },
  {
    label: 'TIDE odds-ratio table',
    href: `${base}/tables/Figure2E_TIDE_response_red_bottom_OR.csv`
  }
]
</script>

<style scoped>
.pancancer-immune {
  color: #102a43;
}

.hero-panel,
.download-panel {
  display: flex;
  justify-content: space-between;
  gap: 28px;
  padding: 28px;
  border: 1px solid #d7e7ef;
  border-radius: 18px;
  background: linear-gradient(135deg, #f5fbfc 0%, #eef6ff 100%);
  box-shadow: 0 12px 30px rgba(15, 58, 85, 0.08);
}

.eyebrow,
.panel-tag {
  margin: 0 0 8px;
  color: #008c95;
  font-weight: 700;
  letter-spacing: 0.02em;
}

h2,
h3 {
  margin: 0;
  color: #073763;
}

h2 {
  font-size: 30px;
}

.hero-text {
  max-width: 860px;
  margin: 12px 0 0;
  color: #536579;
  font-size: 16px;
  line-height: 1.7;
}

.summary-grid {
  display: grid;
  grid-template-columns: repeat(3, 110px);
  gap: 12px;
  align-content: center;
}

.summary-card {
  padding: 16px 12px;
  border-radius: 14px;
  background: #ffffff;
  border: 1px solid #cbe4ee;
  text-align: center;
}

.summary-card strong {
  display: block;
  color: #008c95;
  font-size: 24px;
}

.summary-card span {
  display: block;
  margin-top: 4px;
  color: #52657a;
  font-size: 13px;
}

.figure-grid {
  display: grid;
  grid-template-columns: repeat(2, minmax(0, 1fr));
  gap: 22px;
  margin-top: 24px;
}

.figure-card {
  padding: 20px;
  border: 1px solid #dbe8ef;
  border-radius: 16px;
  background: #ffffff;
  box-shadow: 0 10px 24px rgba(15, 58, 85, 0.07);
}

.figure-card-header {
  display: flex;
  justify-content: space-between;
  gap: 16px;
  align-items: flex-start;
}

.figure-card h3 {
  font-size: 20px;
}

.figure-desc {
  min-height: 46px;
  margin: 10px 0 14px;
  color: #52657a;
  line-height: 1.55;
}

.action-row {
  display: flex;
  gap: 8px;
  flex-wrap: wrap;
  justify-content: flex-end;
}

.action-row a,
.download-list a {
  display: inline-flex;
  align-items: center;
  justify-content: center;
  border-radius: 999px;
  padding: 8px 14px;
  background: #073763;
  color: #ffffff;
  text-decoration: none;
  font-weight: 700;
  font-size: 13px;
}

.action-row a:last-child,
.download-list a {
  background: #008c95;
}

.figure-preview {
  border: 1px solid #e5eef4;
  border-radius: 12px;
  background: #fbfdff;
  overflow: hidden;
}

.figure-preview img {
  display: block;
  width: 100%;
  height: auto;
}

.download-panel {
  margin-top: 24px;
  align-items: center;
}

.download-panel p {
  margin: 10px 0 0;
  color: #52657a;
}

.download-list {
  display: flex;
  flex-wrap: wrap;
  gap: 10px;
  justify-content: flex-end;
  max-width: 560px;
}

@media (max-width: 1100px) {
  .hero-panel,
  .download-panel {
    flex-direction: column;
  }

  .summary-grid,
  .figure-grid {
    grid-template-columns: 1fr;
  }

  .download-list,
  .action-row {
    justify-content: flex-start;
  }
}
</style>
