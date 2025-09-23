import { defineDocumentType, makeSource } from "contentlayer2/source-files"
import rehypeHighlight from "rehype-highlight"
import rehypeSlug from "rehype-slug"
import rehypeAutolinkHeadings from "rehype-autolink-headings"
import remarkGfm from "remark-gfm"

export const Announcement = defineDocumentType(() => ({
  name: "Announcement",
  filePathPattern: `announcements/**/*.mdx`,
  contentType: "mdx",
  fields: {
    title: {
      type: "string",
      required: true,
    },
    publishedAt: {
      type: "date",
      required: true,
    },
    summary: {
      type: "string",
      required: true,
    },
    tags: {
      type: "list",
      of: { type: "enum", options: ["活动", "招新", "比赛", "讲座", "通知"] },
      required: true,
    },
    priority: {
      type: "number",
      required: false,
      default: 0,
    },
    expiresAt: {
      type: "date",
      required: false,
    },
    attachments: {
      type: "list",
      of: {
        type: "object",
        fields: {
          name: { type: "string", required: true },
          url: { type: "string", required: true },
        },
      },
      required: false,
    },
  },
  computedFields: {
    url: {
      type: "string",
      resolve: (doc) => `/announcements/${doc._raw.flattenedPath.replace("announcements/", "")}`,
    },
    slug: {
      type: "string",
      resolve: (doc) => doc._raw.flattenedPath.replace("announcements/", ""),
    },
    isExpired: {
      type: "boolean",
      resolve: (doc) => (doc.expiresAt ? new Date(doc.expiresAt) < new Date() : false),
    },
  },
}))

export const Post = defineDocumentType(() => ({
  name: "Post",
  filePathPattern: `posts/**/*.mdx`,
  contentType: "mdx",
  fields: {
    title: {
      type: "string",
      required: true,
    },
    publishedAt: {
      type: "date",
      required: true,
    },
    updatedAt: {
      type: "date",
      required: false,
    },
    excerpt: {
      type: "string",
      required: true,
    },
    tags: {
      type: "list",
      of: { type: "string" },
      required: true,
    },
    authors: {
      type: "list",
      of: { type: "string" },
      required: true,
    },
    coverImage: {
      type: "string",
      required: false,
    },
  },
  computedFields: {
    url: {
      type: "string",
      resolve: (doc) => `/posts/${doc._raw.flattenedPath.replace("posts/", "")}`,
    },
    slug: {
      type: "string",
      resolve: (doc) => doc._raw.flattenedPath.replace("posts/", ""),
    },
  },
}))

export const LearnResource = defineDocumentType(() => ({
  name: "LearnResource",
  filePathPattern: `learn/resources/**/*.mdx`,
  contentType: "mdx",
  fields: {
    title: {
      type: "string",
      required: true,
    },
    updatedAt: {
      type: "date",
      required: true,
    },
    type: {
      type: "enum",
      options: ["文章", "视频", "项目", "数据集", "幻灯片"],
      required: true,
    },
    level: {
      type: "enum",
      options: ["入门", "进阶", "实战"],
      required: true,
    },
    track: {
      type: "string",
      required: true,
    },
    duration: {
      type: "string",
      required: false,
    },
    prerequisites: {
      type: "list",
      of: { type: "string" },
      required: false,
    },
    externalUrl: {
      type: "string",
      required: false,
    },
  },
  computedFields: {
    url: {
      type: "string",
      resolve: (doc) => `/learn/resources/${doc._raw.flattenedPath.replace("learn/resources/", "")}`,
    },
    slug: {
      type: "string",
      resolve: (doc) => doc._raw.flattenedPath.replace("learn/resources/", ""),
    },
  },
}))

export default makeSource({
  contentDirPath: "./content",
  documentTypes: [Announcement, Post, LearnResource],
  mdx: {
    remarkPlugins: [remarkGfm],
    rehypePlugins: [
      rehypeHighlight,
      rehypeSlug,
      [
        rehypeAutolinkHeadings,
        {
          properties: {
            className: ["subheading-anchor"],
            ariaLabel: "Link to section",
          },
        },
      ],
    ],
  },
})
