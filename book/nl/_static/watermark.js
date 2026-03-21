// Per-page watermark text
// Update this map when document status changes
const watermarks = {
  '01_speciale_relativiteit': { text: '"RELEASE CANDIDATE 2"', size: '4rem' },
  // All other pages fall through to the default "DRAFT" from CSS
};

(function () {
  const path = window.location.pathname;
  for (const [page, config] of Object.entries(watermarks)) {
    if (path.includes(page)) {
      document.documentElement.style.setProperty('--watermark-text', config.text);
      if (config.size) {
        document.documentElement.style.setProperty('--watermark-size', config.size);
      }
      break;
    }
  }
})();
