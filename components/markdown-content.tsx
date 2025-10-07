'use client'

import { useEffect } from 'react'

interface MarkdownContentProps {
  content: string
}

export function MarkdownContent({ content }: MarkdownContentProps) {
  useEffect(() => {
    // 为代码块添加复制按钮和语言标识
    const enhanceCodeBlocks = () => {
      const codeBlocks = document.querySelectorAll('pre code')

      codeBlocks.forEach((codeBlock) => {
        const pre = codeBlock.parentElement
        if (!pre) return

        // 添加语言标识
        const language = codeBlock.className.replace('language-', '').replace('hljs', '').trim()
        if (language && !pre.getAttribute('data-language')) {
          pre.setAttribute('data-language', language.toUpperCase())
        }

        // 添加复制按钮
        if (!pre.querySelector('.copy-button')) {
          const button = document.createElement('button')
          button.className = 'copy-button absolute top-2 left-2 bg-gray-800 text-white px-2 py-1 rounded text-sm hover:bg-gray-700 z-10'
          button.textContent = '复制'
          button.onclick = async () => {
            try {
              await navigator.clipboard.writeText(codeBlock.textContent || '')
              button.textContent = '已复制'
              setTimeout(() => {
                button.textContent = '复制'
              }, 2000)
            } catch (err) {
              console.error('复制失败:', err)
            }
          }

          pre.style.position = 'relative'
          pre.appendChild(button)
        }
      })
    }

    enhanceCodeBlocks()
  }, [content])

  return (
    <div
      className="prose prose-gray dark:prose-invert max-w-none"
      dangerouslySetInnerHTML={{ __html: content }}
    />
  )
}