const puppeteer = require('puppeteer');
const fs = require('fs');

(async () => {
  const browser = await puppeteer.launch();
  const page = await browser.newPage();

  const files = fs.readdirSync('./frames').filter(f => f.endsWith('.html'));
  for (const file of files.sort()) {
    await page.goto(`file://${__dirname}/frames/${file}`);
    await page.screenshot({ path: `frames/${file.replace('.html', '.png')}` });
  }

  await browser.close();
})();
