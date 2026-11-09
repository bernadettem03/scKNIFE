// Prefer light mode on first use of this documentation theme; retain later choices.
(function () {
  try {
    if (!localStorage.getItem("scknife-furo-theme-initialized")) {
      if (!localStorage.getItem("theme")) localStorage.setItem("theme", "light");
      document.body.dataset.theme = localStorage.getItem("theme");
      localStorage.setItem("scknife-furo-theme-initialized", "true");
    }
  } catch (_) { /* Storage may be disabled in private browser settings. */ }
})();
