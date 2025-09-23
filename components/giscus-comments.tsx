"use client"

import { useEffect, useRef } from "react"
import { useTheme } from "next-themes"

interface GiscusCommentsProps {
  repo: string
  repoId: string
  category: string
  categoryId: string
  mapping?: "pathname" | "url" | "title" | "og:title" | "specific" | "number"
  term?: string
  reactionsEnabled?: boolean
  emitMetadata?: boolean
  inputPosition?: "top" | "bottom"
  lang?: string
  loading?: "lazy" | "eager"
}

export function GiscusComments({
  repo = "xiongsircool/sbc-website", // 替换为你的 GitHub 仓库
  repoId = "R_kgDONJQqVw", // 替换为你的仓库 ID
  category = "General", // 替换为你的讨论分类
  categoryId = "DIC_kwDONJQqV84CkQHZ", // 替换为你的分类 ID
  mapping = "pathname",
  term,
  reactionsEnabled = true,
  emitMetadata = false,
  inputPosition = "bottom",
  lang = "zh-CN",
  loading = "lazy",
}: GiscusCommentsProps) {
  const ref = useRef<HTMLDivElement>(null)
  const { theme, resolvedTheme } = useTheme()

  // 根据当前主题确定 Giscus 主题
  const giscusTheme = resolvedTheme === "dark" ? "dark" : "light"

  useEffect(() => {
    if (!ref.current || ref.current.hasChildNodes()) return

    const scriptElem = document.createElement("script")
    scriptElem.src = "https://giscus.app/client.js"
    scriptElem.async = true
    scriptElem.crossOrigin = "anonymous"

    scriptElem.setAttribute("data-repo", repo)
    scriptElem.setAttribute("data-repo-id", repoId)
    scriptElem.setAttribute("data-category", category)
    scriptElem.setAttribute("data-category-id", categoryId)
    scriptElem.setAttribute("data-mapping", mapping)
    if (term) scriptElem.setAttribute("data-term", term)
    scriptElem.setAttribute("data-strict", "0")
    scriptElem.setAttribute("data-reactions-enabled", reactionsEnabled ? "1" : "0")
    scriptElem.setAttribute("data-emit-metadata", emitMetadata ? "1" : "0")
    scriptElem.setAttribute("data-input-position", inputPosition)
    scriptElem.setAttribute("data-theme", giscusTheme)
    scriptElem.setAttribute("data-lang", lang)
    scriptElem.setAttribute("data-loading", loading)

    ref.current.appendChild(scriptElem)
  }, [
    repo,
    repoId,
    category,
    categoryId,
    mapping,
    term,
    reactionsEnabled,
    emitMetadata,
    inputPosition,
    lang,
    loading,
    giscusTheme,
  ])

  // 当主题变化时更新 Giscus 主题
  useEffect(() => {
    const iframe = document.querySelector<HTMLIFrameElement>("iframe.giscus-frame")
    if (!iframe) return

    iframe.contentWindow?.postMessage(
      {
        giscus: {
          setConfig: {
            theme: giscusTheme,
          },
        },
      },
      "https://giscus.app",
    )
  }, [giscusTheme])

  return (
    <div className="mt-8">
      <div className="border-t border-border pt-8">
        <h3 className="text-lg font-semibold mb-4">评论讨论</h3>
        <div ref={ref} />
      </div>
    </div>
  )
}
