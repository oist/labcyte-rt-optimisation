-- https://stackoverflow.com/a/49396058/5309786

-- links-to-markdown.lua
function Link(el)
  el.target = string.gsub(el.target, "%.md", ".html")
  return el
end
