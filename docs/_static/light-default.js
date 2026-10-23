// Prefer light mode on first use of this documentation theme; retain later choices.
(function () {
  try {
    if (!localStorage.getItem("scknife-furo-theme-initialized")) {
      localStorage.setItem("theme", "light");
      localStorage.setItem("scknife-furo-theme-initialized", "true");
    }
  } catch (_) { /* Storage may be disabled in private browser settings. */ }
})();
