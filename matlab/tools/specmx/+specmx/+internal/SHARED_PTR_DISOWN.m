function v = SHARED_PTR_DISOWN()
  persistent vInitialized;
  if isempty(vInitialized)
    vInitialized = specmx(0, 0);
  end
  v = vInitialized;
end
