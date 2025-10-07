import { marked } from 'marked'
import hljs from 'highlight.js'
import fs from 'fs'
import path from 'path'

// 配置 marked
marked.setOptions({
  highlight: function(code, lang) {
    if (lang && hljs.getLanguage(lang)) {
      try {
        return hljs.highlight(code, { language: lang }).value
      } catch (err) {
        console.error('Highlight.js error:', err)
      }
    }
    return hljs.highlightAuto(code).value
  },
  breaks: true,
  gfm: true
})

export async function renderMarkdownToHTML(content: string): Promise<string> {
  try {
    // 使用 marked 渲染 markdown
    let html = marked.parse(content)

    // 后处理：转换 bilibili 链接为 iframe
    let processedHTML = html

    // 处理 bilibili 视频链接格式: https://www.bilibili.com/video/BVxxxx
    processedHTML = processedHTML.replace(
      /https:\/\/www\.bilibili\.com\/video\/(BV\w+)(?:\?t=(\d+(?:\.\d+)?))?/g,
      (match: string, bvid: string, _startTime: string | undefined) => {
        const embedUrl = `//player.bilibili.com/player.html?bvid=${bvid}&page=1&high_quality=1&danmaku=0`
        return `
          <div class="relative aspect-video w-full rounded-lg overflow-hidden bg-black my-6">
            <iframe
              src="${embedUrl}"
              class="w-full h-full"
              scrolling="no"
              border="0"
              frameBorder="no"
              framespacing="0"
              allowFullScreen
            ></iframe>
          </div>
        `
      }
    )

    // 处理已有的 iframe 标签，确保正确的样式
    processedHTML = processedHTML.replace(
      /<iframe src="\/\/player\.bilibili\.com\/player\.html\?bvid=(BV\w+)"[^>]*><\/iframe>/g,
      (match: string, bvid: string) => {
        return `
          <div class="relative aspect-video w-full rounded-lg overflow-hidden bg-black my-6">
            <iframe
              src="//player.bilibili.com/player.html?bvid=${bvid}"
              class="w-full h-full"
              scrolling="no"
              border="0"
              frameBorder="no"
              framespacing="0"
              allowFullScreen
            ></iframe>
          </div>
        `
      }
    )

    return processedHTML
  } catch (error) {
    console.error('Markdown 渲染错误:', error)
    return `<div class="error">内容渲染失败</div>`
  }
}

// 解析frontmatter的简单函数
function parseFrontmatter(content: string) {
  const frontmatterRegex = /^---\s*\n([\s\S]*?)\n---\s*\n([\s\S]*)$/
  const match = content.match(frontmatterRegex)

  if (!match) {
    return {
      frontmatter: {},
      content: content
    }
  }

  const frontmatterStr = match[1]
  const markdownContent = match[2]

  const frontmatter: any = {}

  // 简单的frontmatter解析
  frontmatterStr.split('\n').forEach(line => {
    const colonIndex = line.indexOf(':')
    if (colonIndex > 0) {
      const key = line.substring(0, colonIndex).trim()
      let value = line.substring(colonIndex + 1).trim()

      // 移除引号
      if (value.startsWith('"') && value.endsWith('"')) {
        value = value.slice(1, -1)
      }

      // 处理数组
      if (value.startsWith('[') && value.endsWith(']')) {
        value = value.slice(1, -1)
          .split(',')
          .map(item => item.trim().replace(/^"|"$/g, ''))
          .filter(item => item.length > 0)
      }

      frontmatter[key] = value
    }
  })

  return {
    frontmatter,
    content: markdownContent
  }
}

export async function loadMarkdownContent(contentPath: string) {
  try {
    const fullPath = path.join(process.cwd(), contentPath)

    if (!fs.existsSync(fullPath)) {
      return {
        frontmatter: {},
        content: '<p>内容未找到</p>',
        raw: ''
      }
    }

    const rawContent = fs.readFileSync(fullPath, 'utf-8')
    const { frontmatter, content } = parseFrontmatter(rawContent)

    // Convert Markdown content to HTML
    const html = await renderMarkdownToHTML(content)

    return {
      frontmatter,
      content: html,
      raw: content
    }
  } catch (error) {
    console.error('加载内容错误:', error)
    return {
      frontmatter: {},
      content: '<p>加载失败</p>',
      raw: ''
    }
  }
}

// 为保持兼容性，保留原有函数名
export async function renderMDXToHTML(content: string): Promise<string> {
  return renderMarkdownToHTML(content)
}

export async function loadPostContent(slug: string) {
  const contentPath = `content/posts/${slug}.md`
  const result = await loadMarkdownContent(contentPath)
  return {
    title: result.frontmatter.title || slug,
    content: result.content,
    raw: result.raw
  }
}

export async function loadAnnouncementContent(slug: string) {
  const contentPath = `content/announcements/${slug}.md`
  const result = await loadMarkdownContent(contentPath)
  return {
    title: result.frontmatter.title || slug,
    content: result.content,
    raw: result.raw
  }
}

export async function loadLearnResourceContent(slug: string) {
  const contentPath = `content/learn/resources/${slug}.md`
  const result = await loadMarkdownContent(contentPath)
  return {
    title: result.frontmatter.title || slug,
    content: result.content,
    raw: result.raw
  }
}

export async function loadMeetingContent(slug: string) {
  const contentPath = `content/meetings/${slug}.md`
  const result = await loadMarkdownContent(contentPath)
  return {
    title: result.frontmatter.title || slug,
    content: result.content,
    raw: result.raw
  }
}